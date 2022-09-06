#! /usr/bin/env python3

# Heavily relies on the Illuminaio code https://github.com/HenrikBengtsson/illuminaio/blob/develop/R/readIDAT_nonenc.R


from __future__ import print_function
from   typing import Set
import argparse
import os
import sys
import struct
from  numpy import empty, uint32,fromfile,uint16
import gzip
import re
# we avoid the use of backslashes to assist in templatising the code for Nextflow
TAB=chr(9)
EOL=chr(10)

FID_nSNPsRead    =   1000   
FID_IlluminaID   =    102   
FID_SD           =    103   
FID_Mean         =    104   
FID_NBeads       =    107   
FID_MidBlock     =    200   
FID_RunInfo      =    300   
FID_RedGreen     =    400   
FID_MostlyNull   =    401   
FID_Barcode      =    402   
FID_ChipType     =    403   

FIDs = [FID_nSNPsRead, FID_IlluminaID, FID_SD, FID_Mean, FID_NBeads, FID_MidBlock, FID_RunInfo, FID_RedGreen, FID_MostlyNull, FID_Barcode, FID_ChipType]

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('sample', type=str, metavar='samplesheet'),
    parser.add_argument('idat', type=str, metavar='IDATDIR',help="directory where IDAT files can be found"),
    parser.add_argument('manifest', type=str, metavar='MANIFESTFILE',help="file with Illumina manifest"),
    parser.add_argument('--no-sample-sheet',dest="nosamples",default=False,\
                        action="store_true",\
                        help="use if no samplesheet")
    parser.add_argument('--sample-sample-column',dest="sample_col",type=int,default=0,\
                        help="col in sample file where sample ID can be found (number cols from 0)")
    parser.add_argument('--sex-column',dest="sex_col",type=int,default=0,\
                        help="col in sample file where sex can be found (number cols from 0)")
    parser.add_argument("--sentrix-barcode-column",dest="barcode_col",type=int,default=3,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sentrix-position-column",dest="position_col",type=int,default=4,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sample-delimiter",dest="sample_delimiter",type=str,default=",",\
                        help="what separates entries in the sample file"),
    parser.add_argument("-a","--allow-missing-idat",dest="allow",action="store_true",default=False,\
                        help="if IDAT files are missing, report only, otherwise crash")
    parser.add_argument("-s","--suppress-warnings",dest="suppress",action="store_true",default=False,\
                        help="suppress warnings -- be careful")
    parser.add_argument("--skip-header",action="store_true",dest="has_header",default=False)
    parser.add_argument("-n","--num-threads",dest="num_threads",type=int,default=1,\
                        help="number of threads for parallelism")
    parser.add_argument("-c","--chrom-pos",dest="chrom_pos",action="store_true",default=False,\
                        help="show snp as chromosome-position (default SNP ID as in manifest")
    parser.add_argument("-o","--out",dest="out",type=str,required=True,\
                        help="name of output file")
    args = parser.parse_args()
    return args


class SampleEntry:
    
    def __init__(self,pid,fs):
        self.pid=pid    # person or sample ID
        self.fs = fs    # red and gree files



class SNP:

    def __init__(self,addr_a,addr_b,strand,name,uid,chrom,pos,the_snps):
        self.addr_a   = addr_a
        self.addr_b   = addr_b
        self.strand   = strand
        self.name     = name
        self.uid      = uid
        self.chrom    = chrom
        self.pos      = pos
        self.alleles  = the_snps
        self.a_pos = self.b_pos = None

    def setPos(self,pos,ab): # position of probes in idat file
        if ab == 0:
            self.a_pos = pos
        else:
            self.b_pos = pos


    def chrom_pos(self):
        return (self.chrom,self.pos)

    def __str__(self):
        return("{}:{}".format(self.chrom,self.pos))


def getiDatHash(idat_dirs):
    ''' dict of all idat files and thei locations : different projects organise their 
        idat files differently -- some flat, some in a hierarchical space '''
    hash = {}
    for idat_dir in idat_dirs.split(","):
       tree = os.walk(idat_dir)
       for (d,subds,fs) in tree:
           for f in fs:
               (fn,base)=os.path.splitext(f)
               if f.endswith(".idat") or f.endswith(".idat.gz"):
                   key = fn if base==".gz" else f
                   hash[key] = os.path.join(d,f)
    print (hash)
    return hash


bcwell = {}


def parseSampleSheet(args):
    #parse the sample file to extract the IDs of the particpants and their corresponding
    #idat files. Print warning or crash if the files don't exst
    with open(args.sample) as mf:
        idats=[]

        for line in mf:
            recs    = line.strip().split(args.sample_delimiter)
            pid     = recs[args.sample_col]
            barcode = recs[args.barcode_col]
            pos     = recs[args.position_col]
            sex     = recs[args.sex_col]
            curr_fs = []
            ok= True
            warning = ""
            for colour in ["Red","Grn"]:
                base_file = "{barcode}_{pos}_{colour}.idat".format(barcode=barcode,pos=pos,colour=colour)
                f =  idat_hash[base_file]
                this_ok = os.access(f,os.R_OK)
                if not this_ok: warning=warning+"Warning: file {} does not exist or readable{}".format(f,EOL)
                ok = ok & this_ok
                curr_fs.append(f)
            if not ok:
                if args.allow: 
                    if not args.suppress:
                       sys.stderr(warning+EOL)
                    continue
                else:
                    sys.exit("Missing idat files: "+EOL+warning)
            idats.append(SampleEntry(pid,curr_fs))
            bcwell["{barcode}_{pos}".format(barcode=barcode,pos=pos)]=[pid,sex]
    return idats


def colsOfManifest(fnames):
    ''' return the index(base 0) of the column in the  manifest file for the key fields we need'''
    fields = []
    for name in ["IlmnStrand","Name","SNP","AddressA_ID","AddressB_ID","Chr","MapInfo"]:
        fields.append(fnames.index(name))
    return fields

def getManifest(args):
    # Returns a list of all the SNPs plus an index for each probe saying which SNP it belongs go
    snp_manifest = []
    address_index= {}
    with open(args.manifest) as f:
        line=f.readline()
        while 'IlmnID' not in line:
            line=f.readline()
            print(line)
        fnames = line.split(",")
        cols=colsOfManifest(fnames)
        oldpos=oldchrom=1
        i=0
        for line in f:
           i=i+1
           if i%100000 == 0 : print(i)
           fields   = line.split(",")
           if "Controls" in fields[0]: break
           try:
              (strand,name,snps,address_a,address_b,chrom,pos)=\
                  map(lambda col: fields[col],cols)
              uid = "{}:{}".format(chrom,pos)
              the_snps = snps[1:-1]
              addr_a   = int(address_a)
              addr_b = int(address_b) if address_b else 0
              snp_manifest.append(SNP(addr_a,addr_b,strand,name,uid,chrom,pos,the_snps))
           except IndexError:
              if not args.suppress: sys.stderr.write(line)
        snp_manifest.sort(key=SNP.chrom_pos)
        for i, snp in enumerate(snp_manifest):
            address_index[snp.addr_a]=(i,0)
            if snp.addr_b>0: address_index[snp.addr_b]=(i,1)
    return (snp_manifest,address_index)

def getManifest(args):
    # Returns a list of all the SNPs plus an index for each probe saying which SNP it belongs go
    x=Popen("grep -n Controls] %s | cut -d : -f 1"%args.manifiest,\
            shell=True, stdout=PIPE)
    n=int(x.stdout.readline().decode("utf-8"))-9
    mf = pd.read_csv(args.manifest,skiprows=7,delimiter=",",nrows=n,\
                     dtype={'Chr':str})

def getNum(f,num_bytes=4):
    #j=i+num_bytes
    if num_bytes==2:
        code='H'
    elif num_bytes==4:
        code='L'
    elif num_bytes==8:
        code='Q'
    data = f.read(num_bytes)
    res, = struct.unpack("<%s"%code,data)
    return res

def my_open(f):
    if '.gz' in f:
        return gzip.open(f,'rb')
    else:
        return open(f,'rb')

def getVals(fname):
    
    with my_open(fname) as f:
        # read as string
        magic_number=f.read(4)
        if magic_number.decode('utf-8') != "IDAT":
            print("Magic number is ",magic_number)
            sys.exit("File <%s> not an IDAT file"%fname)
        version = getNum(f,8)
        if  version != 3:
            sys.exit("IDAT version 3 supported only, found {}".format(version))
        fcount = getNum(f)
        field_val = {}
        for i in range(fcount):
            fcode = getNum(f,2)
            offset= getNum(f,8)
            field_val[fcode]=offset
        f.seek(field_val[FID_nSNPsRead])
        num_markers =  getNum(f)
        f.seek(field_val[FID_Barcode])
        bcode       =  getNum(f)
        f.seek(field_val[FID_IlluminaID])
        iids =  fromfile(f,dtype=uint32,count=num_markers)
        f.seek(field_val[FID_Mean])
        vals =  fromfile(f,dtype=uint16,count=num_markers)
        return (iids,vals)


def probeIndexInData(idatf,smf,aidx):
    probe_addr, intensities = getVals(idatf)
    j=0
    for i, addr in enumerate(probe_addr):
        try:
            snp_pos,ab = aidx[addr]
            smf[snp_pos].setPos(i,ab)
            #print(i,addr,snp_pos,ab)
        except KeyError:
            if not args.suppress:
                sys.stderr.write("Warning: Probe {} not in manifest{}".format(addr,EOL))
            j=j+1
    print("There are ",j," beads not in the manifest")
    for snp in smf:
        if not snp.a_pos and not args.suppress:
                sys.stderr.write("Warning: SNP {} not in idat file{}".format(addr,EOL))



def getSNPIntensities(data,s_idx,res,smf):
    ''' For each SNP find probe(s) for that SNP and get values
        data -- numpy array where data to be stored
        sample -- whose file we're dealing with (index)
        res    -- array of red, green values by probe
        smf    -- SNP manififest file '''
    for i, snp in enumerate(smf):
        a_idx = snp.a_pos
        if not a_idx: continue  # warning given earlier
        for colour in [0,1]:
            data[s_idx,i,colour] =  res[colour][a_idx]


        
def showHeading(f,idats):
    f.write("SNP{}Coord{}Alleles".format(TAB,TAB,TAB))
    for entry in idats:
        f.write("{}{}".format(TAB,entry.pid))
    f.write(EOL)

def showSNP(data,snp_i,snp,address_pos,AB):
   if not address_pos:  # no B
       return
   dir1=os.path.join(args.out,snp.chrom)
   the_pos = "%0d"%(int(snp.pos)/100000)
   dir2=os.path.join(dir1,the_pos)
   if not os.path.exists(args.out):
       os.mkdir(args.out)
   if not os.path.exists(dir1):
       os.mkdir(dir1)
   if not os.path.exists(dir2):
       os.mkdir(dir2)
   name = os.path.join(dir2,snp.name)
   if AB : name = name+"_B"
   f = open("%s.csv"%name,"w")
   num=data.shape[0]
   for sample_i in range(0,num):
      for colour in [0,1]:
        f.write("{}{}".format(TAB,data[sample_i,snp_i,colour]))
      f.write(EOL)
   f.close() 


def showIntensities(args,data,smf,idats):
    for snp_i,snp in enumerate(smf):
        showSNP(data,snp_i,snp,snp.addr_a,"")
        showSNP(data,snp_i,snp,snp.addr_b,"B")
        if snp.addr_b:
           showSNP




    
def batchProcessIDATS(args,data,idats,smf):
    g=open("%s_meta.csv"%args.out,"w")
    for i in range(len(idats)):
        sample = idats[i].fs
        m = re.search(".*/(.*)_Red.ida.*",sample[0])
        plate_well = m.group(1)
        if args.sample=='NONE':
            the_id="plate_well\tXX"
        else:
            the_id="%s\n"%"\t".join(bcwell[plate_well])
        g.write(the_id)
        res = list(map(lambda fn : getVals(fn)[1], sample))
        getSNPIntensities(data,i,res,smf)
    g.close()


def processIDATS(args,idats,smf,aidx):
     n_samples=len(idats)
     n_snps   = len(smf)
     data     = empty((n_samples,n_snps,2), dtype=uint32)
     print("Creating an array of size %dx2x%d=%dGB\n"%\
           (n_samples,n_snps,data.nbytes/1024/1024/1024))
     print(idats[1])
     print(idats[1].fs)
     probeIndexInData(idats[1].fs[0],smf,aidx)
     batchProcessIDATS(args,data,idats,smf)
     showIntensities(args,data,smf,idats)



if __name__ == '__main__':
    args       =  parseArguments()
    idat_hash=getiDatHash(args.idat)
    print(list(idat_hash.keys()))
    if args.sample=="NONE":
        bases : Set[str] =set()
        idats=[]
        for k in idat_hash.keys():
            print (k)
            m=re.search("(.*)_(Red|Grn).*", k)
            bases.add(m.group(1))
        for b in bases:
            cfiles = SampleEntry(b,[idat_hash[b+"_Red.idat"], idat_hash[b+"_Grn.idat"]])
            idats.append(cfiles)
    else:
        print("Reading sample sheet")
        idats      =  parseSampleSheet(args) 
    print("Reading manifest")
    (smf,aidx) =  getManifest(args)
    print("Processing idat files")
    processIDATS(args,idats,smf,aidx)
