import requests
import argparse

parser = argparse.ArgumentParser(description='Download results from GNPS workflows')
parser.add_argument('-j', '--job', type=str, help='GNPS job ID')
parser.add_argument('-t', '--type', type=str, help='Workflow type: AreaUnderCurve, SpectralCounts and Network')

args = parser.parse_args()
jobid = args.job

#jobid = '395828bf825f4df78eb8dd8417a69d41'
#jobid = '4443d168a5224a02abe71a13d34a116e'
if args.type=="AreaUnderCurve":
# Feature table, area under the curve
	url = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=" + jobid + "&view=download_consensus_spectra"
elif args.type=="Network":
# Network
	url = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=" + jobid + "&view=download_clustered_spectra"
elif args.type=="SpectralCounts":
# Feature table, spectral counts 
	url = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=" + jobid + "&view=download_cluster_buckettable"

response = requests.post(url)

with open("download.zip", "wb") as output_file:
    for block in response.iter_content(1024):
        output_file.write(block)
