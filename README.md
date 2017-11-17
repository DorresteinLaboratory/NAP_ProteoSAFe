### Network Annotation Propagation

Network Annotation Propagation (NAP) uses spectral networks to propagate information from spectral library matching, in order to improve in silico fragmentation candidate structure ranking.

The ranking improvement is performed using the Fusion concept, proposed on [MetFusion](http://onlinelibrary.wiley.com/doi/10.1002/jms.3123/abstract) publication.

#### Install 

Before installing NAP, install [ProteoSAFe](https://github.com/CCMS-UCSD/ProteoSAFe).

Make sure ProteoSAFe is running properly, you can learn more with the [workflow developer resources](http://proteomics.ucsd.edu/Software/ProteoSAFe/developer_resources/).

Futher documentation can be found at [Workflow XML file documentation](https://bix-lab.ucsd.edu/display/PS/XML+Configuration+Overview).

After ProteoSAFe installation download [MetFrag CL](http://c-ruttkies.github.io/MetFrag/) and place the jar file on the nap_ccms2/Snap directory.
Please refer to [MetFusion](http://onlinelibrary.wiley.com/doi/10.1002/jms.3123/abstract) and [MetFrag](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-148) publications to learn more about the tools.

Edit the file install_all.txt to change the YOURPATHCONDAENV and YOURPATHTOOLSFOLDER variables. Alternatively, refer to ProteoSAFe documentation to encode the path at the tools xml file and avoid hardcoding.

After that, run the installation script to install the conda environment and all dependencies.

```
./install_all.txt
```
After installation move the nap_ccms2 folder to proteosafe/tools folder and nap_ccms2_workflow to proteosafe/workflows.
On proteosafe/workflows rename nap_ccms2_workflow to nap_ccms2.

NAP was installed and tested under 

```
Linux 4.9.0-2-amd64 #1 SMP Debian 4.9.18-1 (2017-03-30) x86_64 GNU/Linux
```

