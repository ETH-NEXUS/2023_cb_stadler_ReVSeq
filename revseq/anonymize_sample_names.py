import pandas as pd
import sys, os
import argparse
import uuid
import base64
import unittest
import gzip

class ethidGenerator(unittest.TestCase):
    @staticmethod
    def encode_base32(uuid4):
        base32uuid = base64.b32encode(uuid4.bytes)
        return base32uuid.decode().replace('=', '')
    @staticmethod
    def decode_base32(ethid):
        return uuid.UUID(bytes=base64.b32decode(ethid + ('=' * (-len(ethid) % 8))))
    def test_eth_id(self):
        uuid4 = uuid.uuid4()
        eth_id = EthIdTestCase.encode_base32(uuid4)
        print(eth_id)
        self.assertEqual(uuid4, EthIdTestCase.decode_base32(eth_id))

def verify_sample(samplename,):
	with gzip.open(sample, 'rb') as f:
    	for i, l in enumerate(f):
        	pass
print("File {1} contain {0} lines".format(i + 1, myfile))

def generate_new_ethids(anon_df, all_samples):
	if len(anon) == 0:
		print("Warning: empty anonymization table. This is just a notice, not an error.")
	newsamples = []
	for samplename in all_samples:
		if samplename not in anon_df["sample_name"]:
			newethid = ethidGenerator.encode_base32(uuid.uuid4())
			while newethid in anon_df["ethid"]:
				newethid = ethidGenerator.encode_base32(uuid.uuid4())
			newsamples.append( [samplename.name, newethid] )
	newsamples = pd.DataFrame(newsamples)
	return newsamples


def link_anonimised_names(samplename, ethid, sampledir, anonymizeddir):
	print("Log: working on ethid: ", ethid)
	sampledir = sampledir + "/" + str(samplename)
	anonymizeddir = anonymizeddir + "/" + str(ethid)
	if not os.path.exists(sampledir):
		sys.exit("Error: there is no directory to be used as origin for the anonymized link for ethid " + str(ethid) + ". This should not happen. Please check if there are problems with the storage.\nPLEASE NOTE that all links for the previous ethids in this run have been already generated and need to be cleaned up before the next run. At this point, the anonymization matching table has not been updated yet")
	try:
		os.mkdir(anonymizeddir)
		for file in os.scandir(sampledir):
			temp_split = file.name.split("_")
			temp_split[0] = str(ethid)
			anonfile = "_".join([temp_split[0], temp_split[2], temp_split[3], temp_split[4]])
			os.symlink(sampledir + "/" + file.name, anonymizeddir + "/" + anonfile)
	except:
		for anonfile in os.scandir(anonymizeddir):
			os.remove(anonymizeddir + "/" + anonfile)
		os.rmdir(anonymizeddir)
		print("Warning: Failed to link the raw data to the anonymized folder for ethid " + str(ethid) + ". The sample will be skipped and not saved in the anonymization table")
		return [samplename, ethid, 1]
	return [samplename, ethid, 0]


def create_sample_map(processed, anondir):
	ethids = [ str(id[1]) for id in processed]
	lanes = []
	for id in ethids:
		rawfiles = os.scandir(anondir + "/" + id)
		lanes_partial = []
		for file in rawfiles:
			lanes_partial.append(file.name.split("_")[1])
		lanes_partial = list(set(lanes_partial))
		lanes_partial = [ [id, lane] for lane in lanes_partial]
		for lane in lanes_partial:
			lanes.append(lane)
	lanes = pd.DataFrame(lanes)
	lanes = lanes.rename(columns={0: "sample", 1: "lane"})
	return lanes


# Script
if __name__ == '__main__':
	# Parse input args
	parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
				  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--anontable', required=True, type=str, help='the table containing the full anonymization list')
	parser.add_argument('--sampledir', required=True, type=str, help='the directory containing all raw data samples, one folder per sample')
	parser.add_argument('--samplemapfile', required=True, type=str, help='the path and filename to the samplemap used by the pipeline')
	parser.add_argument('--anonymizeddir', required=True, type=str, help='the directory containing all anonymized samples, as links to the original data')

	args = parser.parse_args()
    
	anon = pd.read_table(args.anontable, header=0)
	all_samples = os.scandir(args.sampledir)

	newsamples = generate_new_ethids(anon, all_samples)

	linked = [ link_anonimised_names(sample, ethid, args.sampledir, args.anonymizeddir) for sample,ethid in zip(newsamples[0], newsamples[1]) ]
	
	processed = [ sample[0:2] for sample in linked if sample[2]==0 ]
	try:
		pd.DataFrame(processed).to_csv(args.anontable, sep="\t", mode="a", header=None, index=None)
	except:
		sys.exit("Error: could not write the anonymization table. All directories have been created and need to be deleted for the procedure to move forward\nTo know what directories to delete, check in the anonymized directory the ethids that are not included in the anonymised table and delete all of them")
	
	samplemap = create_sample_map(processed, args.anonymizeddir)
	try:
		samplemap.to_csv(args.samplemapfile, sep="\t", mode="a", header=None, index=None)
	except:
		sys.exit("Error: could not write the samplemap table. All directories and new lines in the anonymization table have been created and need to be deleted for the procedure to move forward.")
	
	
