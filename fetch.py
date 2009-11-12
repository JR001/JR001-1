#!/usr/local/bin/python
from Bio import Entrez
from Bio import SeqIO
import Bio.Align.Applications  #for dealing with Muscle command line calls 

##Function to take list of accessions and download the data from Entrez
# Function to take list of accessions and download the data from Entrez
# Function to convert that .gb file to a custom FASTA file.
# Modified from: http://www.biopython.org/wiki/BioSQL
# and Julius Lucks' script

#Notes on usage.  To run it, use the two functions at the very bottom.
#If you are downloading new files from Genbank, then use fetch.download, give it a list of accession numbers in a text file, one line per numbers
#If you are parsing a Genbank file and creating a FASTA file, use fetch.parse_gb_to_FASTA and give it
#your accession list and the output file name. (you must have already downloaded these using download, and run it in the same dir)
#however, the lists do not have to be the same, or in the same order

############# set your file name and dir here, this is the only one you need to change, rest are based on that.
filename_base = "/Users/johncumbers/Desktop/41_Chro_gb/Chro_accessions_3prime_intact"  # name and loc of your input accesssion number list this.txt
#############
# no need to touch these below:
input_file_name_list_of_accession_numbers = "%s.txt" % filename_base 
output_file_name_FASTA = "output_%s.FASTA" % filename_base
output_file_name_FASTA_aligned =  "output_%s.FASTA.aligned" % filename_base


#
# Things to do: add some feedback to this instruction, print out as each one is downloaded
def download(accession_list):  #download each from genbank, save each as a .gb file in current dir
	try:
		for line in open(accession_list):  # opened in text-mode; all EOLs are converted to '\n'
			try:
				if not line.startswith("#"):  #see below for details of partition function and what it does to deal with #commenting
					line = line.partition(' #')[0]  # comments start with a space#
					line = line.rstrip('\n')
					handle = Entrez.efetch(db="nuccore", id=line, rettype="gb")
					genbank_string = handle.read()
					gb_file_name = line+'.gb'
					f = open(gb_file_name, 'w')
					f.write(genbank_string)
					print "Created file named %s" % gb_file_name
	#################				f.write('//\n')  #not sure why/if we want to include this Genbank separator but Julius did
					f.close()
			except IOError:
				print 'Network error, are you connected to the internet?'
				quit()
				raise
	except IOError:
		print 'something wrong with the input file? Please check the spelling and location'
		quit()
		raise


#WARNING: if there is a wrong accession number then the program will not tell you about it, and just download a blank file
#You'll only know after you try and parse the file with the parse function.
def parse_gb_to_FASTA(accession_list, fasta_filename):   #parse each of the files you just  
#downloaded and outputs them in a customized FASTA file format.
	output_handle = open(fasta_filename, "w")  #setup the FASTA output file
	try:
		parsed = [ ] #this could be used to return everything that was passed
		for line in open(accession_list):	#read each line one by one of the input file, loop through them until end of file
			if not line.startswith("#"):
				line = line.partition(' #')[0]  # Comments start with a space#
				
				#this effectively ignores anything after a hash on each line, 
				#allowing comments to be added to the text file at the beginning of lines or after accesssion numbers
				# anything after a hash will be ignore and the next line will be processed
				#partition returns a tuple: everything before the partition string, the partition string, 
				#and everything after the partition string. So, by indexing with [0] we take just the part before the partition string.
				
				line = line.rstrip('\n')			#removed all the end of line characters
				gb_file_name = line + '.gb'			 # we previously downloaded accession list as .gb files
				
				try:
					gb_file = file(gb_file_name, 'r')  #read each file
				except IOError:
					print 'Is the file %s downloaded?' % gb_file_name
					raise
				gb_parsed_record = SeqIO.parse(gb_file, "genbank").next()  
				# seqIO method converts the genbank file .gb into a generator
				# we can use the next() command on the generator to retrieve 
				# an object representing the parsed file as above.
				gb_file.close()
				# now we get to do stuff with the data inside the .gb file
				# so let's write out the bits that we want.
				output_handle.write(">%s %s\n%s\n" % (  # > var var end line var end line. (approx. FASTA format)
					gb_parsed_record.id, 
					gb_parsed_record.description, 
					gb_parsed_record.seq.tostring().lower())) 
			print "%s" % line	#Print all lines and comments in the file, so help you see what's going on
			
		#print 'FASTA file saved to:\n%s\n' % fasta_filename 
	except IOError:
		print 'Is the file %s downloaded and have you checked what is inside it?  Maybe tiy got the accesssion number wrong? I am quitting' % line
		quit()
		raise


	
''' #NOT USED
handle = Entrez.efetch(db="nuccore", id="6273291,6273290,6273289", rettype="gb")
for seq_record in SeqIO.parse(handle, "gb") :
    print seq_record.id, seq_record.description[:50] + "..."
    print "Sequence length %i," % len(seq_record.seq),
    print "from: %s" % seq_record.annotations['source']
handle.close()
'''

''' #NOT USED
for seq_record in SeqIO.parse(handle, "gb") :
    print seq_record.id, seq_record.description[:50] + "..."
    print "Sequence length %i," % len(seq_record.seq),
    print "from: %s" % seq_record.annotations['source']
handle.close()
'''

#  Use these parts here just add your file names in place.
#  run in TextMate etc.  with command + R.
#  Comment or uncomment the parts you want to run.
#
#import fetch  # no need to do this if you are running it from within.
#fetch.download('/Users/johncumbers/Desktop/41_Chro_gb/41Chro_accessions.txt')
#fetch.download('/Users/johncumbers/Desktop/41_Chro_gb/outgroup_sequences.txt')


#fetch.parse_gb_to_FASTA('/Users/johncumbers/Desktop/41_Chro_gb/41Chro_accessions.txt', '/Users/johncumbers/Desktop/41_Chro_gb/41Chro_accessions.FASTA')
#fetch.
#parse_gb_to_FASTA('/Users/johncumbers/Desktop/41_Chro_gb/Chro_16s_ends_intact_with_outgroups.txt', '/Users/johncumbers/Desktop/41_Chro_gb/Chro_16s_ends_intact_with_outgroups.FASTA')

#parse_gb_to_FASTA('input_file_name_list_of_accession_numbers', '/Users/johncumbers/Desktop/41_Chro_gb/Chro_accessions_3prime_intact.FASTA')
parse_gb_to_FASTA('/Users/johncumbers/Desktop/41_Chro_gb/Chro_accessions_5prime_intact.txt', 'output_file_name_FASTA')


