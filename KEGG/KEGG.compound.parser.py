import Bio.KEGG.Compound

KEGGw = open('KEGGyTheKEGG.fq','w')

handle = open("compound")
for record in Bio.KEGG.Compound.parse(handle):
	#print record.name[0]
	#print record.formula
	KEGGw.write("%s\t%s\n" % (record.name[0], record.formula))