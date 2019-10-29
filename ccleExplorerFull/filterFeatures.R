#a = load('/mnt/efs/data/BigF_19q3_20190820_interpolated_cors.RData')
#b = load('/mnt/efs/data/BigF_19q3_20190820_zeroSDremoved.RData')
rm(list = ls()); gc(); 
a = load('/mnt/efs/data/BigF_19q3_20190928_interpolated_cors.RData')
b = load('/mnt/efs/data/BigF_19q3_20190928_zeroSDremoved.RData')



foi = as.integer(Fo)

#tfoi = table(foi)
tfoi = rep(0, ncol(Fo));
N = length(foi)
cat(N,'\n')
for(i in 1:N){
  if(i%%1000000==0){cat(i,' ')}
  j = foi[i]; 
  tfoi[j]=tfoi[j]+1; 
}

#ttfoi = table(table(foi))
ttfoi = table(tfoi)

sort(ttfoi)[1:20]
-sort(-ttfoi)[1:20]

names(tfoi)= colnames(F)

sort(-tfoi)[1:10]


###



flset = c('exonUsuge','meth','sanger.CN','expr','sanger.expr','DEMETER2','fusionUnfilt','mut','geneloss','fusion','mut.aa')

foi2 = as.integer(Fo[1:50,which(!(F.type %in% flset))])

#tfoi = table(foi)
tfoi2 = rep(0, ncol(Fo));
N = length(foi2)
cat(N,'\n')
for(i in 1:N){
  if(i%%100000==0){cat(i,' ')}
  j = foi2[i]; 
  tfoi2[j]=tfoi2[j]+1; 
}

ttfoi2 = table(tfoi2)

names(tfoi2)= colnames(F)

sort(-tfoi2)[1:50]

sort(-ttfoi2)[1:10]

a
b

table(F.type[(tfoi2==0)|(tfoi==1)])
ifilt = which(!((tfoi2==0)|(tfoi==1)))

F = F[,ifilt]
F.type = F.type[ifilt]
F.gene = F.gene[ifilt]
Fo = Fo[,ifilt]
Fc = Fc[,ifilt]


#save(F,F.type,F.gene, depmap.id,  file = '/mnt/efs/data/BigF_19q3_20190820_zeroSDremoved_flt300k.RData')
#save(Fo,Fc, file = '/mnt/efs/data/BigF_19q3_20190820_interpolated_cors_flt300k.RData')
save(F,F.type,F.gene, depmap.id,  file = '/mnt/efs/data/BigF_19q3_20190928_zeroSDremoved_flt300k.RData')
save(Fo,Fc, file = '/mnt/efs/data/BigF_19q3_20190928_interpolated_cors_flt300k.RData')

#load( file = '/mnt/efs/data/BigF_19q3_20190820_zeroSDremoved_flt300k.RData')
#load( file = '/mnt/efs/data/BigF_19q3_20190820_interpolated_cors_flt300k.RData')
Fof = Fo;
Fcf = Fc;
fnam = colnames(F)

#load( file = '/mnt/efs/data/BigF_19q3_20190820_interpolated_cors.RData')
#load( file = '/mnt/efs/data/BigF_19q3_20190820_zeroSDremoved.RData')
load( file = '/mnt/efs/data/BigF_19q3_20190928_interpolated_cors.RData')
load( file = '/mnt/efs/data/BigF_19q3_20190928_zeroSDremoved.RData')
jnew = match(colnames(F), fnam)

Fofn = matrix(jnew[Fof], nrow = nrow(Fof))
Fo = Fofn;
Fc = Fcf;
#save(Fo,Fc, file = '/mnt/efs/data/BigF_19q3_20190820_interpolated_cors_flt300k.RData')
save(Fo,Fc, file = '/mnt/efs/data/BigF_19q3_20190928_interpolated_cors_flt300k.RData')
