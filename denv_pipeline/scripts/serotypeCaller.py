#!/usr/bin/env python

import os
import sys

from Bio import SeqIO

def main():

    minCov=50

    if os.path.exists(sys.argv[1]):
        tmp=[str(i).strip() for i in open(sys.argv[1],"r")]

        if(len(tmp)>1):

            for i in sys.argv[1:]:
                a=SeqIO.parse(i,'fasta')

                for i,k in enumerate(a):
                    if i==0:
                        seqLen=len(k.seq)
                        seqLenNoAmb=len(str(k.seq).replace('n','').replace('-','').replace('N',''))
                        percCov=round(seqLenNoAmb/seqLen*100,2)

                        if len(sys.argv[1:])==2:
                            if os.path.exists(sys.argv[2]):
                                trimPos=[int(r) for r in [str(j).strip().split("\t") for j in open(sys.argv[2],"r")][0]]
                                seqTrim=str(k.seq[int(trimPos[0])-1:int(trimPos[1])])
                                seqLenTrim=len(seqTrim)
                                seqLenNoAmbTrim=len(seqTrim.replace('n','').replace('-','').replace('N',''))
                                percCovTrim=round(seqLenNoAmbTrim/seqLenTrim*100,2)
                                
                                if percCovTrim>=minCov:
                                    print(str(k.id)+'\t'+str(seqLen)+'\t'+str(seqLenNoAmb)+'\t'+str(percCov)+'\t'+str(percCovTrim))
                                else:
                                    print( "/".join(str(k.id).split("/")[0:2]) +'\t'+'NA'+'\t'+'NA'+'\t'+str(seqLen)+'\t'+str(seqLenNoAmb)+'\t'+str(percCov)+'\t'+str(percCovTrim))

                            trimpFile=open( str(sys.argv[1]).replace(".out.aln",".out.trim.aln"),"w")
                            trimpFile.write(">"+str(sys.argv[1]).split(".")[0]+"\n")
                            trimpFile.write(seqTrim+"\n")
                            trimpFile.close()

                        else:
                            if percCov>=minCov:
                                print(str(k.id)+'\t'+str(seqLen)+'\t'+str(seqLenNoAmb)+'\t'+str(percCov)+'\t'+str('NA'))
                            else:
                                print("/".join(str(k.id).split("/")[0:2]) +'\t'+'NA'+'\t'+'NA'+'\t'+str(seqLen)+'\t'+str(seqLenNoAmb)+'\t'+str(percCov)+'\t'+str('NA'))

        else:
            print(str(sys.argv[1]).split(".")[0]+'/'+str(".".join(str(sys.argv[1]).split(".")[0:3])+".cons.fa" )+'/Unknown\t'+'NA'+'\t'+'NA'+'\t'+'NA'+'\t'+str(0)+'\t'+str(0))
    else:
        print(str(sys.argv[1]).split(".")[0]+'/'+str(".".join(str(sys.argv[1]).split(".")[0:3])+".cons.fa" )+'/Unknown\t'+'NA'+'\t'+'NA'+'\t'+'NA'+'\t'+str(0)+'\t'+str(0))

if __name__=="__main__":
    main()


