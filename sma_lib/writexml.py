#take in parameter object (defined in parameters.py) and output xml file.
#into file named in outxml
#import xml.etree.cElementTree as ET
from yattag import indent
from yattag import Doc
import datetime
# based on http://stackoverflow.com/questions/3605680/creating-a-simple-xml-file-using-python	
#and http://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
#FIXME: unfortunately, this is in alphabetical order instead of nicely organized.
#to get a readable form, use http://www.yattag.org/ 's indent function
# def write_xml0(par,outxml):
	# root = ET.Element('root')
	# doc = ET.SubElement(root,'doc')
	
	#get elements of par object. 	#remove values like __class__
	# elts = [a for a in dir(par) if not a.startswith('__')]

	# for elt in elts:
		# val = getattr(par,elt)
		# t = type(val).__name__
		# ET.SubElement(doc,elt,type=t).text = str(val)
	
	# tree = ET.ElementTree(root)
	# tree.write(outxml)

#with reasonable formatting
def write_xml(par,outxml,antype,file,outcomes):
	doc,tag,text=Doc().tagtext()
	
	#get elements of par object. 	#remove values like __class__
	elts = [a for a in dir(par) if not a.startswith('__')]
	
	#notes
	with tag('analysis notes'):
		with tag("Analysis Type"):
			text(antype)
		with tag("Completed at"):
			text(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
		with tag("on file"):
			text(file)
		#with tag('in folder'):
		#	text(folder)
		
	#outcomes, for ffpdax
	if(antype =="ffpdax"):
		with tag('outcomes'):
			with tag('number of peaks'):
				text(str(outcomes[0]))
			with tag('mean duration'):
				text(str(outcomes[1]))
			with tag('median duration'):
				text(str(outcomes[2]))
		
	#parameters
	with tag("Settings"):
		for elt in elts:
			val = getattr(par,elt)
			t = type(val).__name__
			#print elt,t,val
			with tag(elt+ ' type = \''+t+"'"):
				text(str(val))
				
	formed = indent(doc.getvalue())
	#print 'xml:'
	#print doc.getvalue()
	f = open(outxml,'w')
	f.write(formed)
	