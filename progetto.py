#!/usr/bin/python
# :-)

"""
Codice per l'esame
"""

#P1
"""
dato un dizionario che associa gli id del FASTA file alla stringa di DNA, ritorna il numero di record nel dizionario
"""
def get_num_records(records):
    """codice da sviluppare"""
    return len(records)

#P2
"""
dato un dizionario che associa gli id del FASTA file alla stringa di DNA, ritorna un dizionario che associa gli id del dizionario alla lunghezza della stringa.
"""
def get_length_map(records):
    """codice da sviluppare"""
    return {key:len(value) for key,value in records.items()}


"""
dato un dizionario che associa gli id alla lunghezza della stringa, ritorna la lunghezza della stringa minima e la lunghezza della stringa massima
"""
def get_min_max_length(len_map):
    """codice da sviluppare"""
    return min(len_map.values()),max(len_map.values())

"""
dato un dizionario che associa gli id alla lunghezza della stringa, ritorna il numero delle sequenze che hanno lunghezza minima e massima
"""
def get_min_max_num(len_map):
    """codice da sviluppare"""
    min,max=get_min_max_length(len_map)
    values=list(len_map.values())
    return (values.count(min),values.count(max))

"""
dato un dizionario che associa gli id alla lunghezza della stringa, ritorna la lista di id associati alle sequenze che hanno lunghezza minima e massima
"""
def get_min_max_length_id(len_map):
    """codice da sviluppare"""
    min,max = get_min_max_length(len_map)
    minresult,maxresult = [ key for key,value in len_map.items() if value == min],[key for key,value in len_map.items() if value == max]
    return (minresult,maxresult)


#P3
"""
data una sequenza di dna ed un frame di lettura, restituisce la lista di tutti gli ORF contenuti nella stringa di DNA.
"""
def get_ORF_list(dna,frame=0):
    dna_split = [dna.lower()[x:x+3] for x in range(frame,len(dna)-3+1,3)]
    result=[]
    for i,x in enumerate(dna_split):
        if x == 'atg':
            for j,y in enumerate(dna_split[i:]):
                if y in {'tga','tag','taa'}:
                    result.append(''.join(dna_split[i:i+j+1]))
                    break
            else:
                break

    return result

"""
dato un dizionario che associa gli id del FASTA file alla stringa di DNA, restituisce un dizionario che associa gli id alla lista di ORF contenuti nella sequenza corrispondente.
"""
def get_ORF_map(records,frame=0) :
    """codice da sviluppare"""
    return {key:get_ORF_list(value,frame) for key,value in records.items()}

"""
DAto un dizionario che associa gli id alla lista degli ORF per quella sequenza, restituisce l'id associato all'ORF piu' lungo, l'ORF piu' lungo e la lunghezza di tale ORF (in caso ci fossero piu' ORF massimali restituire uno qualsiasi di questi ORF).
"""
def longest_ORF(orf_map):
    """codice da sviluppare"""
    orf_map_max = { key:max([ (len(x),x) for x in values],default=(0,'')) for key,values in orf_map.items() } # comprimo le liste negli orf piu' lunghi
    max_len,max_orf,max_id = max([(values[0],values[1],key) for key,values in orf_map_max.items() ],default=(0,'',''))
    return max_id,max_orf,max_len

#P4
"""
data una lunghezza n ed una lista di sequenze di nucleotidi, restituisce un dizionario che associa tutte le sottostringhe ripetute di lunghezza n (i.e., sottostringhe che compaiono almeno una volta in almeno una sequenza) al loro numero di occorrenze in tutte le sequenze.
"""
def get_all_repeats_aux(n,word): # returns dictionary for all n-length substrings in word
    rep_list = [word.lower()[x:x+3] for x in range(len(word)-3+1)]
    return {x:rep_list.count(x) for x in set(rep_list)}

def get_all_repeats(n,seq_list):
    dict_list = [get_all_repeats_aux(n,word) for word in seq_list]
    return {x:sum([d[x] for d in dict_list if x in d]) for x in set.union(*[set(d) for d in dict_list])}

"""
dato un dizionario che associa sottostringhe al numero di occorrenze, ritorna la sottostringa con il numero massimo di occorrenze (ed il numero di occorrenze)
"""
def most_frequent_repeat(rep_map):
    """codice da sviluppare"""
    max_occorrenze, max_string = max([ (value,key) for key,value in rep_map.items() ])
    return max_string,max_occorrenze

"""
dato un dizionario che associa sottostringhe al numero di occorrenze ed un numero di occorrenze occ, ritorna l'insieme delle stringhe che ha un numero di occorrenze associato pari almeno a occ.
"""
def filter_repeats(rep_map,occ):
    """codice da sviluppare"""
    return [ key for key,value in rep_map.items() if value >= occ ]


if __name__ == "__main__" :
    filename = input("Scegliere file fasta:\n\ta) dna.simple.fasta\n\tb) dna.long.fasta\n")
    if (filename == 'a'):
        filename = 'dna.simple.fasta'
    else:
        filename = 'dna.long.fasta'

    from fastautil import read_fasta

    import sys

    #test P1
    records = read_fasta(filename)
    print("Test numero record, a) 5 b) 28 ")
    print("numero record = %s " % get_num_records(records))

    #test P2
    len_map = get_length_map(records)
    #print("dizionario lunghezze = %s " % len_map)
    min_seq,max_seq = get_min_max_length(len_map)
    print("test lunghezza a) min: 29 max: 112 b) min: 186 max: 1686")
    print("lunghezza seq minima %s e massima %s " % (min_seq,max_seq))
    n_min,n_max = get_min_max_num(len_map)
    print("numero lunghezza seq a) min: 2 max: 1 b) min: 1 max: 1")
    print("numero seq con lunghezza minima %s e massima %s " % (n_min,n_max))
    id_min,id_max =    get_min_max_length_id(len_map)
    print("id seq a) min: ['id2','id3'] max: ['id4'] b) min: ['lcl|EU819142.1_cds_ACF06958.1_18'] max: ['lcl|EU819142.1_cds_ACF06952.1_12']")
    print("id seq con lunghezza minima %s e massima %s " % (id_min,id_max))


    #test P3
    test_string = 'CATGCTATGGTGCCCTAAAAGATGACGCCCTACCCCCCCCCTAGTGATTAGTTTCGAGACATACTGTGTTT'
    print("test_string contiene 2 ORF in frame 0 ed 1 in frame 1")
    print("lista ORF associata alla stringa %s : %s" % (test_string,get_ORF_list(test_string)))
    orf_map    = get_ORF_map(records)
#    print("dizionario ORF = %s " % orf_map)
    print("numero ORF totali: a) 2 b) 253 ")
#    print("numero ORF totali: %s " % orf_map.values())
    print("numero ORF totali: %d " % sum([len(k) for k in orf_map.values()]))
    id_lorf,lorf,length_lorf = longest_ORF(orf_map)
#    print("id sequenza ORF piu' lungo = %s, ORF piu' lungo = %s, lunghezza ORF = %d" % (id_lorf,lorf,length_lorf))
    print("id ORF piu' lungo  e lunghezza orf a) id4, 54 b) lcl|EU819142.1_cds_ACF06952.1_12, 1686")
    print("id sequenza ORF piu' lungo = %s, lunghezza ORF = %d" % (id_lorf,length_lorf))


    #Test P4
    n = 3
    rep_map = get_all_repeats(n,records.values())
    #print('dizionario delle sottostringhe ripetute di lunghezza %d : %s' % (n,rep_map))

    m_f_rep,rep_n = most_frequent_repeat(rep_map)
    print("sottostringa ripetuta piu' di frequente e numero di occorrenze a) CCC, 70 b) GGC, 954")
    print("sottostringa ripetuta piu' di frequente nel dizionario = %s , numero di occorrenze = %d" % (m_f_rep,rep_n))

    occ = 30
    rep_set = filter_repeats(rep_map,occ)
    print(rep_set)
    print("dimensione insieme a) 1 b) 64 ")
    print("dimensione insieme di tutte le sottostringhe ripetute che occorrono almeno %d volte = %s " % (occ,len(rep_set)))
