Instructions to download the mouse embryo data at time 8.0 days, 8.25, and 8.50 days after fertilization:

First, download the file here https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6967/E-MTAB-6967.sdrf.txt

Then, the links to download the various files can be extracted for the various timepoints using the commands below:


8.0 days post-fertilization:

`awk -F "\t" '{if($1=="Sample 16")  {print $54"\n"$56"\n"$58}}' E-MTAB-6967.sdrf.txt > links.txt`


8.25 days post-fertilization:

`awk -F "\t" '{if($1=="Sample 28")  {print $54"\n"$56"\n"$58}}' E-MTAB-6967.sdrf.txt > links_25.txt`


8.50 days post-fertilization:

`awk -F "\t" '{if($1=="Sample 17")  {print $54"\n"$56"\n"$58}}' E-MTAB-6967.sdrf.txt > links_5.txt`
