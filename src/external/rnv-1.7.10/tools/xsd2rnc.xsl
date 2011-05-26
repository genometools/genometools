<?xml version="1.0"?>

<xsl:transform 
	xmlns:xsd="http://www.w3.org/2001/XMLSchema"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	version="1.0">

  <xsl:output method="text"/>

  <xsl:template match="/">
    default namespace = "<xsl:value-of select="/xsd:schema/xsd:element/@name"/>-NS" 
    namespace xsi = "http://www.w3.org/2001/XMLSchema-instance"

    start = element <xsl:value-of select="/xsd:schema/xsd:element/@name"/> { 
      attribute xsi:schemaLocation { text } ? ,
      <xsl:apply-templates select="/xsd:schema/xsd:simpleType[1]/xsd:restriction"/>
    }
  </xsl:template>  

  <xsl:template match="xsd:restriction">
      xsd:<xsl:value-of select="@base"/> {
	<xsl:for-each select="*">
	   <xsl:value-of select="name()"/>="<xsl:value-of select="@value"/>"
	</xsl:for-each>
      }
  </xsl:template>

</xsl:transform>
