"""
modulo che contiene funzioni di utilita' per file in formato fasta
"""

def read_fasta(filename) :
	"""
	legge tutte le sequenze contenute in un 
	file fasta e associa gli id alle sequenze 
	"""
	seqs={}
	try :
		f = open(filename)
	except IOError:
		print("the file does not exists")
		return seqs

	for line in f :
		#toglie \n se presente
		line = line.rstrip()
		#controlla se header
		if line.startswith('>') :
			words=line.split() #divide una stringa
			name=words[0][1:]
			seqs[name]=''
		else : #se non e' un header
			seqs[name] = seqs[name] + line

	f.close()	
	return seqs	

def print_fasta_seqs(filename):
	"""
	assume filename is a fasta file, reads all sequences in a dictionary and prints the dictionary
	"""
	seqs = read_fasta(filename)
	for name,seq in seqs.items() :
		print(name,seq)


