# $Id: rvp.py,v 1.5 2004/11/09 11:29:31 dvd Exp $

# embedding sample for RVP, a part of RNV, http://davidashen.net/rnv.html
# code kept simple to show the technique, not to provide a general purpose
# module.
#
# details of the protocol are in a long comment near the start of rvp.c
#

import sys, os, string, re, xml.parsers.expat

global rvpin,rvpout,pat,errors,parser,text,ismixed,prevline,prevcol

# raised by resp if it gets something the module does not understand
class ProtocolError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

# run RVP with grammar specified on the command line
def launch():
  global rvpin,rvpout
  inp,out=os.popen2('rvp '+string.join(sys.argv[1:],' '))
  rvpin,rvpout=inp.fileno(),out.fileno() # os.read|os.write will be used

# terminate string with zero, encode in utf-8 and then send to RVP
def send(s):
  os.write(rvpin,s.encode('UTF-8')+'\0')

# receive a zero-terminated response from RVP, zero byte is dropped
def recv():
  s=''
  while 1:
    s=s+os.read(rvpout,16) # 16 is good for ok responses, errors should be rare
    if(s[-1]=='\0'): break # last character in a response is always '\0'
  return s[:-1]

# return current pattern value, and print error message on stderr, if any
def resp():
  global errors,prevline,prevcol
  r=string.split(recv(),' ',3)
  if(r[0]=='ok'): return r[1]
  if(r[0]=='error'):
    errors=1
    if(r[3]!=''): # if the error string is empty, then error
    		  # is occured in erroneous state; don't report
      line,col=parser.ErrorLineNumber,parser.ErrorColumnNumber
      if(line!=prevline or col!=prevcol): # one report per file position
	sys.stderr.write(str(line)+","+str(col)+": "+r[3])
	prevline,prevcol=line,col
    return r[1]
  if(r[0]=='er'):
    errors=1
    return r[1]
  raise ProtocolError,"unexpected response '"+r[0]+"'"

def start_tag_open(cur,name):
  send('start-tag-open '+cur+' '+name)
  return resp()

def attribute(cur,name,val):
  send('attribute '+cur+' '+name+' '+val)
  return resp()

def start_tag_close(cur,name):
  send('start-tag-close '+cur+' '+name)
  return resp()

def end_tag(cur,name):
  send('end-tag '+cur+' '+name)
  return resp()

def textonly(cur,text):
  send('text '+cur+' '+text)
  return resp()

# in mixed content, whitespace is simply discarded, and any
# non-whitespace is equal; but this optimization gives only
# 5% increase in speed at most in practical cases
def mixed(cur,text):
  if(re.search('[^\t\n ]',text)):
    send('mixed '+cur+' .')
    return resp()
  else: return cur

def start(g):
  send('start '+g)
  return resp()

def quit():
  send('quit')
  return resp()

# Expat handlers

# If I remember correctly, Expat has the custom not to concatenate
# text nodes; therefore CharDataHandler just concatenates them into
# text, and then flush_text passes the text to the validator
def flush_text():
  global ismixed,pat,text

  if(ismixed):
    pat=mixed(pat,text)
  else:
    pat=textonly(pat,text)
  text=''

def start_element(name,attrs):
  global ismixed,pat

  ismixed=1
  flush_text()
  pat=start_tag_open(pat,name)
  ismixed=0
  for n,v in attrs.items(): pat=attribute(pat,n,v)
  pat=start_tag_close(pat,name)

def end_element(name):
  global ismixed,pat

  flush_text()
  pat=end_tag(pat,name)
  ismixed=1

def characters(data):
  global text

  text=text+data

# Main

errors=0
launch()
pat=start('0') # that is, start of the first grammar;
	       # multiple grammars can be passed to rvp
parser=xml.parsers.expat.ParserCreate('UTF-8',':') # last colon in the name
                                  # separates local name from namespace URI
parser.StartElementHandler=start_element
parser.EndElementHandler=end_element
parser.CharacterDataHandler=characters

text=''
prevline,prevcol=-1,-1
parser.ParseFile(sys.stdin)

quit() # this stops RVP; many files can be validated with a single RVP
       # running, concurrently or sequentially
sys.exit(errors)
