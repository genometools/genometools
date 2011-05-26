# $Id: gn.sed 275 2004-01-09 20:38:01Z dvd $
# convert head of arx to arguments of rvp
/grammars/,/}/ {
  s/[ 	]+/ /g
  s/#.*$|(^| )grammars( |{|$)|{|}/ /g
  s/ *([^ ]+)="([^ ]+)" */\1=\2 /gp
}
