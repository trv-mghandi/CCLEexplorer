DEMO_mode=TRUE; ## if DEMO_mode true, it loads a smaller dataset, otherwise loads the whole 600k+ features matrix 


library(shiny)
library(visNetwork)
dataDir = '/home/ubuntu/efs/data'
dataDir = '/mnt/efs/data'
dataDir = '/Users/mghandi/Documents/efs_local_copy/data'
dataDir = '../../ccleExplorerData'

  #(genes, rsub.idx1, rsub.idx2, rsub.pv, file = paste0(ddir,'single/all/edges.Rdata'))
  ddir= paste0(dataDir,'/codep/')
  #load(file = paste0(ddir,'single/all/edges.Rdata'))
  
  nmax.options = c(1:29,10*(3:20));  #max number of neighbors
  dmax.options = c(1);#,2)
  emax.options = c(1:30,40,50,100, 500,1000); #max number of edges shown 
  
  #edges_netformat = data.frame(id = 1:length(rsub.idx1),from = rsub.idx1, to = rsub.idx2, value = -log10(rsub.pv), color =    adjustcolor( '#97C2FC', alpha.f = 0.8))
  #nodes_netformat = data.frame(id = 1:length(genes), title = genes, label = genes)
  ####

  enets <- new.env()
  enets$network_idx = array(dim = c(length(nmax.options),length(dmax.options),length(emax.options)))   ## this is the index for edges and nodes list objects for different values of Nmax, Dmax(max distance), emax(maximum edges to show)
  enets$edges_list = list(); 
  enets$nodes_list = list(); 
  enets$listSize=0; 
  enets$geneName = ""
  enets$featureName = ""
  enets$featureIdx = 1
  

  
  if(FALSE){
    library(broman)
    
    f.types.cnt = table(F.type)
    f.types.id = names(f.types.cnt)
    f.types.nam = paste0(f.types.id, ' (',f.types.cnt,')')
    #f.types.id
    f.types.cols = matrix(c("annot", 'Laser Lemon',
                            "Avana", 'Sky Blue',#Blue Green',
                            "DEMETER2", 'Indigo',
                            'Sanger.Score',"Robin's Egg Blue",
                            'Sanger.drug','Fern',#Asparagus',
                            'CCLE.drug','Magic Mint',
                            'PRISM.drug',"Screamin' Green",
                            'CTDD.drug','Jungle Green',
                            'exonUsuge','Desert Sand',
                            'expr','Orange Red',
                            'tpm','Scarlet',
                            'sanger.expr','Jazzberry Jam',
                            'miRNA','Razzle Dazzle Rose',
                            'ssGSEA','Salmon',
                            'sanger.pathway','Lavender',
                            'CN','Violet (Purple)',
                            'sanger.CN','Vivid Violet',
                            'meth','Neon Carrot',
                            'meth.enh','Sunglow',
                            'mut','Brown',
                            'mut.aa','Shadow',#Outrageous Orange',
                            'fusion','Sepia',
                            'fusionUnfilt','Chestnut',
                            'geneloss','Outer Space',
                            'metabolites','Peach',#Eggplant',
                            'chrom','Gold',#Beaver',
                            'RPPA','Wisteria'), ncol =2, byrow = TRUE)
    
    f.types.cols[which(is.na(match(f.types.cols[,2],names( brocolors("crayons"))))),]
    f.types.cols[,2] =brocolors("crayons")[f.types.cols[,2]]
    colnames(f.types.cols) = c('f.types.id','col' )
    
    
    plot(1:nrow(f.types.cols), rep(0, nrow(f.types.cols)), pch=20, cex=3, col=f.types.cols[,2])             
    text((1:nrow(f.types.cols))-0.5, rep(0, nrow(f.types.cols)),paste(' ',f.types.cols[,1]), srt=90, pos = 4)          
    
    save(    f.types.cnt, 
             f.types.id,
             f.types.nam,
             f.types.cols, 
             file = '/mnt/efs/data/f.types.cols.RDAT'
    )
  }
  #load('/mnt/efs/data/f.types.cols.RDAT')
  load(paste0(dataDir,'/f.types.cols.RDAT'))
  
  f.type.selected = f.types.nam[grep('AVANA', toupper(f.types.nam))[1]];
  f.type.id.selected = f.types.id[match(f.type.selected, f.types.nam)]
  
  g.cur.f.types.nam = f.types.id;
  g.cur.f.types.nam.selected = g.cur.f.types.nam ; #all selected by default 
  g.cur.f.types.nam.selected = setdiff(g.cur.f.types.nam.selected, 'exonUsuge')
  
  
#  f.fnams.in.type = colnames(F)[F.type == f.type.id.selected]
  
  
    
  if(DEMO_mode){
    #a=load( file ='/mnt/efs/data/BigF_sample_test.RDAT')
    load(paste0(dataDir,'/DEMO/BigF_sample_test.RDAT'))
    
  }else{
  #  load(paste0(dataDir,'/BigF_19q3_20190820_interpolated_cors.RData'))
  #  load(paste0(dataDir,'/BigF_19q3_20190820_zeroSDremoved.RData'))
    load(paste0(dataDir,'/Full/BigF_19q3_20190928_interpolated_cors_flt300k.RData'))
    load(paste0(dataDir,'/Full/BigF_19q3_20190928_zeroSDremoved_flt300k.RData'))
  }
  
  calcNetworkBigF = function(featureName='Avana:BRAF', nmax=50, dmax=1, emax=100, f.type.ids.filt =f.types.id, enets){  ## return network id; results in enets$edges_list and enets$nodes_list
    #browser()
    
    enets$curIdx = 1; 
    
    if(nmax>nrow(Fo)){
      nmax = nrow(Fo)
    }

    if(enets$featureName !=featureName){
      j = match(featureName, colnames(F))
      if(is.na(j)){
        return(1)
      } 
      
      enets$network_idx = array(dim = c(length(nmax.options),length(dmax.options),length(emax.options)))   ## this is the index for edges and nodes list objects for different values of Nmax, Dmax(max distance), emax(maximum edges to show)
      enets$edges_list = list(); 
      enets$nodes_list = list(); 
      enets$listSize=0; 
      enets$featureName = featureName
      enets$featureIdx = j
    }
    
    inmax = match(nmax, nmax.options)
    idmax = match(dmax, dmax.options)
    iemax = match(emax, emax.options)
    
    #browser()
    if(!is.na(enets$network_idx[inmax,idmax,iemax])){
      # enets$curIdx = enets$network_idx[inmax,idmax,iemax]
      #  return(enets$network_idx[inmax,idmax,iemax])
      
    }
    
    if(FALSE){
      featureName='Avana:BRAF';
      nmax=10;
      dmax=1; 
      emax=100;
      f.type.selected = f.types.nam[grep('AVANA', toupper(f.types.nam))[1]];
      f.type.id.selected = f.types.id[match(f.type.selected, f.types.nam)]
      
      ii =calcNetworkBigF(featureName, nmax, dmax, emax, f.types.id, enets)
    }  

    
    enets$listSize=  enets$listSize+1; 
    
    
    enets$listSize=1; ######## for now, no hashing!
    
    idxnew =   enets$listSize; #length(enets$edges_list)+1
    enets$edges_list[idxnew]=c()
    enets$nodes_list[idxnew]=c()
    enets$network_idx[inmax,idmax,iemax]=idxnew; 
  #  igi = match(geneName, genes)

    # if(FALSE){    
    #   ## check if we have the solution for a higher V-max
    #   oiemax = which(!is.na(enets$network_idx[inmax,idmax,])); oiemax=oiemax[length(oiemax)]
    #   if(oiemax>iemax){
    #     ## there exists some solution
    #     inetold = enets$network_idx[inmax,idmax,oiemax]
    #     edges_new = enets$edges_list[[inetold]]
    #     nodes_new = enets$nodes_list[[inetold]]
    #     if(length(edges_new)>emax){
    #       ## find list of genes in the top emax connections, and only keeps those (keeps all the edges connecting those; the number could be higher than emax)
    #       nodes_new = c(igi, as.numeric(rbind(rsub.idx1[edges_new[1:emax]],rsub.idx2[edges_new[1:emax]]))); 
    #       nodes_new = nodes_new[!duplicated(nodes_new)]
    #       edges_new = edges_new[(rsub.idx1[edges_new]%in%nodes_new)&(rsub.idx2[edges_new]%in%nodes_new)]
    #     }
    #     enets$edges_list[[idxnew]] <- edges_new
    #     enets$nodes_list[[idxnew]] <- nodes_new
    #     return(idxnew)
    #   }
    #   
    #   ## check if we have the solution for a higher n-max
    #   oinmax = which(!is.na(enets$network_idx[,idmax,iemax])); oinmax=oinmax[length(oinmax)]
    #   if(oinmax>inmax){
    #     ## there exists some solution
    #     inetold = enets$network_idx[oinmax,idmax,iemax]
    #     edges_new = enets$edges_list[[inetold]]
    #     nodes_new = enets$nodes_list[[inetold]]
    #     if(sum(edges_new>nmax)>0){
    #       edges_new = edges_new[edges_new<=nmax]; 
    #       
    #       if(dmax==1){
    #         edges_new_1 = edges_new[(rsub.idx1[edges_new]==igi)|(rsub.idx2[edges_new]==igi)]; ## one side connected 
    #       }else{
    #         inetold = calcNetwork(geneName, nmax, dmax=1, emax,enets)
    #         net1gene_i = enets$nodes_list[[inetold]]
    #         edges_new_1 = edges_new[(rsub.idx1[edges_new]%in%net1gene_i)&(rsub.idx2[edges_new]%in%net1gene_i)]; ## one side connected 
    #       }
    #       
    #       nodes_new = c(igi, as.numeric(rbind(rsub.idx1[edges_new_1],rsub.idx2[edges_new_1]))); 
    #       nodes_new = nodes_new[!duplicated(nodes_new)]
    #       edges_new = edges_new[(rsub.idx1[edges_new]%in%nodes_new)&(rsub.idx2[edges_new]%in%nodes_new)]
    #       
    #     }
    #     enets$edges_list[[idxnew]] <- edges_new
    #     enets$nodes_list[[idxnew]] <- nodes_new
    #     return(idxnew)
    #   }
    # }
    # 
    # 
    
    if(idmax==1){
#      node_ids = Fo[1:nmax, enets$featureIdx]
      node_ids = Fo[, enets$featureIdx]
      
      node_ids = node_ids[F.type[node_ids] %in%f.type.ids.filt]
      if(!(enets$featureIdx %in% node_ids)){
        node_ids = c(enets$featureIdx, node_ids)
      }
      if(length(node_ids)>nmax){node_ids=node_ids[1:nmax]}

      #if(length(node_ids)==1){
      #  browser()
      #}
            
      node_labels = colnames(F)[node_ids]
      node_color = f.types.cols[match(F.type[node_ids],f.types.cols[,1]),2]
      ## need to add the type filter 
      
      fosub = Fo[,node_ids,drop=FALSE]; 
      x.idx = matrix(node_ids, nrow = nrow(fosub), ncol = ncol(fosub), byrow = TRUE)
      
      edges = data.frame(idx1 = as.integer(x.idx), idx2= as.integer(fosub), r = as.numeric(Fc[,node_ids,drop=FALSE]))
      edges = subset(edges, idx2 %in% node_ids)

      
      jj = which(edges$idx2 < edges$idx1)
      if(length(jj)>0){
        hh = edges$idx2[jj]; 
        edges$idx2[jj]=edges$idx1[jj]; 
        edges$idx1[jj]=hh; 
      }
      edges = subset(edges, !duplicated(paste(idx1,idx2,sep=':')))
      edges = subset(edges, idx1!=idx2)
#browser();
      if(nrow(edges)>emax){
        #edges = edges[unique(c(order(-abs(edges$r))[1:emax], which((edges$idx1 == enets$featureIdx)|(edges$idx2 == enets$featureIdx)))),]
        edges = edges[unique(
          c(which((edges$idx1 == enets$featureIdx)|(edges$idx2 == enets$featureIdx)), order(-abs(edges$r))))[1:emax],]
      }
      
      #just for nicer visualization purposes: 
      jj = which(edges$idx2 == enets$featureIdx)
      if(length(jj)>0){
        edges$idx2[jj]=edges$idx1[jj]; 
        edges$idx1[jj]=enets$featureIdx; 
      }

    }
    
    nodes_new = data.frame(id = node_ids, title = node_labels, label = node_labels, color = node_color)
    
    if(nrow(edges)>0){ 
      edges$color = 'steelblue'; edges$color[which(edges$r<0)]='red'
      edges_new = data.frame(id = 1:nrow(edges), from=as.integer(edges$idx1), to=as.integer(edges$idx2), value=abs(edges$r), color=edges$color)
    } else {
      edges_new =subset(data.frame(id=1,from=1,to=1,value=1,color='red'), id==2); 
    }
    
    enets$edges_list[[idxnew]] <- edges_new
    enets$nodes_list[[idxnew]] <- nodes_new
    
    
    enets$curIdx = idxnew
    return(idxnew)
    
  }
  
  
  
  calcNetwork = function(geneName, nmax, dmax, emax, enets){  ## return network id; results in enets$edges_list and enets$nodes_list
    #browser()
    if(FALSE){
      geneName = 'ESR1'
      nmax = 100000;
      dmax = 1;#2; 
      emax = 100; 
      ii =calcNetwork(geneName, nmax, dmax, emax, enets)
    }  
    
    if(enets$geneName !=geneName){
      enets$network_idx = array(dim = c(length(nmax.options),length(dmax.options),length(emax.options)))   ## this is the index for edges and nodes list objects for different values of Nmax, Dmax(max distance), emax(maximum edges to show)
      enets$edges_list = list(); 
      enets$nodes_list = list(); 
      enets$listSize=0; 
      enets$geneName = geneName
    }
    
    
    inmax = match(nmax, nmax.options)
    idmax = match(dmax, dmax.options)
    iemax = match(emax, emax.options)
    
    if(!is.na(enets$network_idx[inmax,idmax,iemax])){
      return(enets$network_idx[inmax,idmax,iemax])
    }
    
    enets$listSize=  enets$listSize+1; 
    idxnew =   enets$listSize; #length(enets$edges_list)+1
    enets$edges_list[idxnew]=c()
    enets$nodes_list[idxnew]=c()
    enets$network_idx[inmax,idmax,iemax]=idxnew; 
    igi = match(geneName, genes)
    
    ## check if we have the solution for a higher V-max
    oiemax = which(!is.na(enets$network_idx[inmax,idmax,])); oiemax=oiemax[length(oiemax)]
    if(oiemax>iemax){
      ## there exists some solution
      inetold = enets$network_idx[inmax,idmax,oiemax]
      edges_new = enets$edges_list[[inetold]]
      nodes_new = enets$nodes_list[[inetold]]
      if(length(edges_new)>emax){
        ## find list of genes in the top emax connections, and only keeps those (keeps all the edges connecting those; the number could be higher than emax)
        nodes_new = c(igi, as.numeric(rbind(rsub.idx1[edges_new[1:emax]],rsub.idx2[edges_new[1:emax]]))); 
        nodes_new = nodes_new[!duplicated(nodes_new)]
        edges_new = edges_new[(rsub.idx1[edges_new]%in%nodes_new)&(rsub.idx2[edges_new]%in%nodes_new)]
      }
      enets$edges_list[[idxnew]] <- edges_new
      enets$nodes_list[[idxnew]] <- nodes_new
      return(idxnew)
    }
    
    ## check if we have the solution for a higher n-max
    oinmax = which(!is.na(enets$network_idx[,idmax,iemax])); oinmax=oinmax[length(oinmax)]
    if(oinmax>inmax){
      ## there exists some solution
      inetold = enets$network_idx[oinmax,idmax,iemax]
      edges_new = enets$edges_list[[inetold]]
      nodes_new = enets$nodes_list[[inetold]]
      if(sum(edges_new>nmax)>0){
        edges_new = edges_new[edges_new<=nmax]; 
        
        if(dmax==1){
          edges_new_1 = edges_new[(rsub.idx1[edges_new]==igi)|(rsub.idx2[edges_new]==igi)]; ## one side connected 
        }else{
          inetold = calcNetwork(geneName, nmax, dmax=1, emax,enets)
          net1gene_i = enets$nodes_list[[inetold]]
          edges_new_1 = edges_new[(rsub.idx1[edges_new]%in%net1gene_i)&(rsub.idx2[edges_new]%in%net1gene_i)]; ## one side connected 
        }
        
        nodes_new = c(igi, as.numeric(rbind(rsub.idx1[edges_new_1],rsub.idx2[edges_new_1]))); 
        nodes_new = nodes_new[!duplicated(nodes_new)]
        edges_new = edges_new[(rsub.idx1[edges_new]%in%nodes_new)&(rsub.idx2[edges_new]%in%nodes_new)]
        
      }
      enets$edges_list[[idxnew]] <- edges_new
      enets$nodes_list[[idxnew]] <- nodes_new
      return(idxnew)
    }
    
    if(idmax==1){
      iisel = which((igi== rsub.idx1[1:nmax])|(igi== rsub.idx2[1:nmax]))
      nodes_new <- as.numeric(rbind(rsub.idx1[iisel],rsub.idx2[iisel])); 
      nodes_new = c(igi, nodes_new[nodes_new!=igi]) 
      edges_new <- which((rsub.idx1[1:nmax] %in%  nodes_new) & (rsub.idx2[1:nmax] %in%  nodes_new))
    }else{
      ## if max dist is 2, then first calc,  maxdist = 1
      if(TRUE){  ## recursive is preferred and is original design, but for some R issues it seems not work in shiny
        inetold = calcNetwork(geneName, nmax, dmax=1, emax,enets)
        net1gene_i = enets$nodes_list[[inetold]]
      }else{  ## non-rec
        iisel = which((igi== rsub.idx1[1:nmax])|(igi== rsub.idx2[1:nmax]))
        nodes_new <- as.numeric(rbind(rsub.idx1[iisel],rsub.idx2[iisel])); 
        net1gene_i = c(igi, nodes_new[nodes_new!=igi]) 
        
      }
      
      iisel = which((rsub.idx1[1:nmax] %in% net1gene_i )|(rsub.idx2[1:nmax] %in% net1gene_i))
      nodes_new <- as.numeric(rbind(rsub.idx1[iisel],rsub.idx2[iisel]));
      nodes_new = nodes_new[!duplicated(nodes_new)]
      edges_new <- which((rsub.idx1[1:nmax] %in%  nodes_new) & (rsub.idx2[1:nmax] %in%  nodes_new))
    }
    
    if(length(edges_new)>emax){
     ## find list of genes in the top emax connections, and only keeps those (keeps all the edges connecting those; the number could be higher than emax)
      nodes_new = c(igi, as.numeric(rbind(rsub.idx1[edges_new[1:emax]],rsub.idx2[edges_new[1:emax]]))); 
      nodes_new = nodes_new[!duplicated(nodes_new)]
      edges_new = edges_new[(rsub.idx1[edges_new]%in%nodes_new)&(rsub.idx2[edges_new]%in%nodes_new)]
    }
    
    enets$edges_list[[idxnew]] <- edges_new
    enets$nodes_list[[idxnew]] <- nodes_new
    return(idxnew)
    
  }
  
  
  
  if(FALSE){
    testfun = function(u,v){
      calcNetwork(geneName = 'FOXA1',10000,1,500,enets)
    }
  
    testfun()
    enets$nodes_list
  }
  
  
  
  
  
  
  
  
  #######
  gi = 'KRAS'
  gi = 'C10orf111'
  gi = 'ESR1'
  
  
     
  library(shiny)
    
    if(FALSE){
      library(devtools)
      library(withr)
    #  with_libpaths(new = "/usr/local/lib/R/site-library", install_github('hadley/shinySignals'))
     # with_libpaths(new = "/usr/local/lib/R/site-library", install_github('jcheng5/bubbles'))
      with_libpaths(install_github('hadley/shinySignals'))
      with_libpaths(install_github('jcheng5/bubbles'))
      install.packages("shinydashboard")
      install.packages(c("shinydashboard",'dplyr','plotly','DT'))
      
    }
    
    library(shinySignals)   # devtools::install_github("hadley/shinySignals")
      
      library(dplyr)
  library(shinydashboard)
  #library(bubbles)        # devtools::install_github("jcheng5/bubbles")

      
  library(plotly)
  library(DT)

      
  
  
  #load(paste0(dataDir,'/Annotations_19q1.RData'))
  load(paste0(dataDir,'/Annotations_19q3.RData'))
  Annotations$inferred_ethnicity = Annotations$genetic_race3
  Annotations$Primary.Disease = Annotations$disease
  Annotations$Subtype.Disease = Annotations$disease_subtype
  Annotations = Annotations[match(rownames(F), rownames(Annotations)),]
  
  
  col_type_table=read.delim(paste0(dataDir,'/col_type_table_19q3.txt'), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  #col_type_tableOld=read.delim(paste0(dataDir,'/col_type_table.txt'), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  
  g = list(); 
  g$Data = list(); 
  
  g$color_type_family = unique(col_type_table$family)
  g$col_type_table =col_type_table
    familySelected = 'type_refined'#'PATHOLOGIST_ANNOTATION'#'type_refined'
      g$col=Annotations[,familySelected] 
      coltab = subset(g$col_type_table, col_type_table$family==familySelected)
  
      g$cols = coltab$col
      
      
      g$cols[which(is.na(g$cols))]="#AAAAAA"
      g$col = factor(g$col, levels = c(coltab$type,'N/A'))
      g$col[which(is.na(g$col))]='N/A'
      
      
      
  g$doubleClick = 0; #to avoid two clicks bug 