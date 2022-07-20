# Tools
## Steps to install GATK
1. Download `gatk-4.x.x.zip` file from [this link](https://github.com/broadinstitute/gatk/releases)
```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.2.3.0/gatk-4.2.3.0.zip
```
2. Unzip the package on remote server
```bash
unzip gatk-4.2.3.zip
```

3. Set up PATH for jar files 
Note: change `/path/to/gatk-package/` to the absolute path of your unzip folder from step 2.
```bash
export PATH="/path/to/gatk-package/":$PATH
```
You can also add this to your `.bashrc` file, so you don't need to run it everytime you want to use GATK. 
```bash
echo "export PATH='/path/to/gatk-package/':$PATH >> ~/.bashrc
source ~/.bashrc
```

4. Test that it works
```bash
gatk info
```
5. Use GATK 
creates necessary reference files for next steps
```bash
gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict
``` 
## Steps to download SnpEff

1. Download latest version
```bash
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
```
2. Unzip file
```bash
unzip snpEff_latest_core.zip
```
3. Configure
```bash
export PATH=$PATH:/path/to/snpEff/
```
4. Download Reference Database
```bash
 java -jar ./snpEff/snpEff.jar download GRCh38.99
```
