###Construct a triplet of lnc or mRNA
ceRNA_mir_ceRNA<-function(ceRNA,ceRNA_list,mir_list){
  #' @param		ceRNA:		candidate mRNA/lncRNA
  #' @param		ceRNA_list:		The interactions of mRNA/lncRNA-miRNA 
  #' @param		mir_list:		The interactions of miRNA-mRNA/lncRNA
  mirs<-ceRNA_list[[ceRNA]]
  common_mirs<-intersect(mirs,names(mir_list))
  temp_list<-mir_list[common_mirs]
  lens<-unlist(lapply(temp_list,length))
  mir_ceRNA<-cbind(rep(names(temp_list),lens),unlist(temp_list))
  l<-dim(mir_ceRNA)[1]
  ceRNA_mir_ceRNAs<-cbind(rep(ceRNA,l),mir_ceRNA)
  rownames(ceRNA_mir_ceRNAs)<-paste(ceRNA_mir_ceRNAs[,1],ceRNA_mir_ceRNAs[,2],ceRNA_mir_ceRNAs[,3],sep = "_")
  
  return(ceRNA_mir_ceRNAs)
}

##Dividing expression profiles into groups based on the state of a genetic alteration,case_normal_profile calls profile_class
##
profile_class<-function(flag,gene_exp,mir_M){
  ##Extracting expression data from three levels(miRNA,lncRNA,mRNA) of a class of flag(cnv state) samples
  ##The result is a list, with the gene/lnc expression profile and the miRNA expression profile.
  #' @param		flag:		a group of sample 
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  a_list<-list()
  length(a_list)<-2
  names(a_list)<-c("gene","mir")
  
  a_list$gene<-gene_exp[,flag]
  a_list$mir<-mir_M[,flag]
  return(a_list)
  
  
}
################################################################
##
case_normal_profile<-function(flag_v,gene_exp,mir_M){
  #' @param		flag_v:		cnv state of a candidate gene
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  #' @param		min_thr:		a sample number threshold. If the number of samples in a certain state is less than min_thr, subsequent operations will be delete.
  
  alter_sample<-tapply(flag_v,as.factor(flag_v),names)##Extract the corresponding sample for each state
  
  result<-lapply(alter_sample,profile_class,gene_exp=gene_exp,mir_M=mir_M)##Extract expression profile
  return(result)
  
}

##apply(cnv_M[1:2,],1,case_normal_profile,gene_M=gene_case,lnc_M=lnc_case,mir_M=mir_case,min_thr=5)
#######################################################################################################################



## Calculate correlation coefficient,corrs calls cors
#######################################################################################
##
cors<-function(a_pair,expr1,expr2){
  #' @param		a_pair:		a vector with 2 RNA names 
  #' @param		expr1和expr2:		Representing the corresponding expression profile
  #print(a_pair)
  temp<-cor.test(expr1[a_pair[1],],expr2[a_pair[2],])
  pcor<-temp$estimate
  p<-temp$ p.value
  result<-c(pcor,p)
  return(result)
  
}
##
corrs<-function(all_pairs,expr1,expr2){
  #' @param		all_pairs:		all RNA names pairs
  #' @param		expr1和expr2:		Representing the corresponding expression profile
  c_p<-t(apply(all_pairs,1,cors,expr1=expr1,expr2=expr2))
  colnames(c_p)<-c("cor","pvalue")
  c_p<-cbind(all_pairs,as.data.frame(c_p))
  return(c_p)
}
########################################################################################

### Significant screen of RNA pairs
pair_filter<-function(c_p,flag,method,cor_thr,pthr){
  #' @param		c_p:		Is a matrix of 4 columns,first two columns are RNA pairs,the 3 column is the coefficient,the 4 column is the p value.
  #' @param		flag:		Negative correlation if it is miRNA
  if(flag=="mir"){##Correlation coefficient screen
    pos1<-which(c_p[,3]<=cor_thr)
  }else{
    pos1<-which(c_p[,3]>=cor_thr)
  }
  
  Padj<-p.adjust(c_p[,4],method=method)
  c_p<-cbind(c_p,Padj)
  Padj<-c_p[,5]
  pos2<-which(Padj<=pthr)##Padj screen
  pos<-intersect(pos1,pos2)
  return(c_p[pos,])
}
########################################################################

##Construct significant ceRNA triplets
filter1<-function(ceRNA,sig_pair,flag){
  ###screen ceRNA
  #' @param		ceRNA:		All candidate triplets
  #' @param		sig_pair:		significant RNA pairs
  #' @param		flag:		which column to screen 
  if(length(ceRNA)==0|length(sig_pair)==0){
    return(NULL)
  }else if(length(ceRNA)==3){
    ceRNA<-matrix(ceRNA,nrow=1)
    
  }
  sig_pairs<-rownames(sig_pair)
  if(flag=="12"){
    temp_pair<-paste(ceRNA[,1],ceRNA[,2],sep="_")
    
  }else if(flag=="13"){
    temp_pair<-paste(ceRNA[,1],ceRNA[,3],sep="_")
  }else{
    temp_pair<-paste(ceRNA[,2],ceRNA[,3],sep="_")
  }
  pos<-which(temp_pair%in%sig_pairs)
  
  return(ceRNA[pos,])
  
}

#############################################################

###Correlation coefficient and p-values of ceRNA
extract_cp<-function(ceRNA,flag,c_ps){
  ###
  #' @param		ceRNA:		All candidate triplets
  #' @param		flag:		which column to screen 
  #' @param		c_ps:		result of pair_filter
  if(length(ceRNA)==0){
    return(NULL)
  }else if(length(ceRNA)==3){
    ceRNA<-matrix(ceRNA,nrow=1)
    
  }
  
  
  
  sig_pairs<-rownames(c_ps)
  
  if(flag=="12"){
    name<-paste(ceRNA[,1],ceRNA[,2],sep="_")
  }else if(flag=="13"){
    name<-paste(ceRNA[,1],ceRNA[,3],sep="_")
  }else{
    name<-paste(ceRNA[,2],ceRNA[,3],sep="_")
  }
  #result<-cbind(ceRNA,c_ps[name,c("cor","pvalue","Padj")])
  result<-cbind(ceRNA,c_ps[name,c("cor","pvalue")])
  return(result)
}
#################################################################

cor_all<-function(ceRNA_triplets,mir_M,gene_exp,cthr1,cthr2,cthr3,pthr){
  #######Produce each RNA pairs
  #' @param		ceRNA_triplets:		result of ceRNA_mir_ceRNA
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		cthr1:		Correlation coefficient threshold,mir与ceRNA1
  #' @param		cthr2:		Correlation coefficient threshold,mir与ceRNA2
  #' @param		cthr3:		Correlation coefficient threshold,ceRNA2与ceRNA1
  #' @param		pthr:		Correlation coefficient significance threshold
  mir_ceRNA1<-unique(ceRNA_triplets[,1:2])
  colnames(mir_ceRNA1)<-c("ceRNA1","mir")
  mir_ceRNA2<-unique(ceRNA_triplets[,c(2,3)])
  colnames(mir_ceRNA2)<-c("mir","ceRNA2")
  ceRNA1_ceRNA2<-unique(ceRNA_triplets[,c(1,3)])
  colnames(ceRNA1_ceRNA2)<-c("ceRNA1","ceRNA2")
  
  ##Calculation correlation coefficient  
  mir_ceRNA1_c_p<-corrs(mir_ceRNA1,expr1=gene_exp,expr2=mir_M)
  mir_ceRNA2_c_p<-corrs(mir_ceRNA2,expr1=mir_M,expr2=gene_exp)
  ceRNA1_ceRNA2_c_p<-corrs(ceRNA1_ceRNA2,expr1=gene_exp,expr2=gene_exp)
  save(mir_ceRNA1_c_p,mir_ceRNA2_c_p,ceRNA1_ceRNA2_c_p,file="c_p.Rdata")
  cmc_ceRNA1_ceRNA2_c_p<-cbind(ceRNA_triplets[rownames(ceRNA1_ceRNA2),],ceRNA1_ceRNA2_c_p[,-(1:2)])
  rownames(cmc_ceRNA1_ceRNA2_c_p)<-paste(cmc_ceRNA1_ceRNA2_c_p[,1],cmc_ceRNA1_ceRNA2_c_p[,2],cmc_ceRNA1_ceRNA2_c_p[,3],sep = "_")
  cmc_ceRNA1_ceRNA2_c_p<-cmc_ceRNA1_ceRNA2_c_p[,-(1:3)]
  save(cmc_ceRNA1_ceRNA2_c_p,file = "cmc_ceRNA1_ceRNA2_c_p.Rdata")
  
  ######Use the correlation filter RNA pairs
  
  mir_ceRNA1_c_p_sig<-pair_filter(mir_ceRNA1_c_p,"mir","fdr",cthr1,pthr)
  rownames(mir_ceRNA1_c_p_sig)<-paste(mir_ceRNA1_c_p_sig[,1],mir_ceRNA1_c_p_sig[,2],sep="_")##rename rownames,use RNA pairs
  
  mir_ceRNA2_c_p_sig<-pair_filter(mir_ceRNA2_c_p,"mir","fdr",cthr2,pthr)
  rownames(mir_ceRNA2_c_p_sig)<-paste(mir_ceRNA2_c_p_sig[,1],mir_ceRNA2_c_p_sig[,2],sep="_")##rename rownames,use RNA pairs
  
  ceRNA1_ceRNA2_c_p_sig<-pair_filter(ceRNA1_ceRNA2_c_p,"gene","fdr",cthr3,pthr)
  rownames(ceRNA1_ceRNA2_c_p_sig)<-paste(ceRNA1_ceRNA2_c_p_sig[,1],ceRNA1_ceRNA2_c_p_sig[,2],sep="_")##rename rownames,use RNA pairs
  
  save(mir_ceRNA1_c_p_sig,mir_ceRNA2_c_p_sig,ceRNA1_ceRNA2_c_p_sig,file="sig_c_p.Rdata")
  #############################################################################################
  #####
  
  miRNA_lncRNA_gene1<-filter1(ceRNA_triplets,mir_ceRNA1_c_p_sig,"12")
  miRNA_lncRNA_gene2<-filter1(miRNA_lncRNA_gene1,mir_ceRNA2_c_p_sig,"23")
  sig_miRNA_lncRNA_gene<-filter1(miRNA_lncRNA_gene2,ceRNA1_ceRNA2_c_p_sig,"13")
  save(sig_miRNA_lncRNA_gene,file="sig_miRNA_lncRNA_gene.Rdata")
  ###############################################################################################
  ##
  
  cp1<-extract_cp(sig_miRNA_lncRNA_gene,"12",mir_ceRNA1_c_p_sig)
  cp2<-extract_cp(cp1,"13",ceRNA1_ceRNA2_c_p_sig)
  sig_miRNA_lncRNA_gene_cp<-extract_cp(cp2,"23",mir_ceRNA2_c_p_sig)
  save(sig_miRNA_lncRNA_gene_cp,file="sig_miRNA_lncRNA_gene_cp.Rdata")
  return(sig_miRNA_lncRNA_gene)
  
}

##
cor_alls<-function(a_flag,exprs,ceRNA_triplets,ceRNA_mir_list,mir_ceRNA_list,cthr1,cthr2,cthr3,pthr){
  #' @param		a_flag:		cnv state
  #' @param		exprs:		result of case_normal_profile
  #' @param		ceRNA_triplets:		result of ceRNA_mir_ceRNA
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		cthr1:		Correlation coefficient threshold,mir与ceRNA1
  #' @param		cthr2:		Correlation coefficient threshold,mir与ceRNA2
  #' @param		cthr3:		Correlation coefficient threshold,ceRNA2与ceRNA1
  #' @param		pthr:		Correlation coefficient significance threshold
  dir.create(a_flag)##create folders for different cnv states of each gene
  temp_path<-getwd();
  temp_path<-paste(temp_path,"/",a_flag,sep="")
  setwd(temp_path)##create a folder for each state of the gene, save the middle results
  
  a_list<-exprs[[a_flag]]
  gene_exp<-a_list[[1]]
  
  mir_M<-a_list[[2]]
  result<-cor_all(ceRNA_triplets,mir_M,gene_exp,cthr1,cthr2,cthr3,pthr)
  
  setwd("../")##
  return(result)
}

ceRNAs_sig<-function(geneid,cnv_M,mir_M,gene_exp,ceRNA_mir_list,mir_ceRNA_list,cthr1,cthr2,cthr3,pthr){
  #' @param		geneid:		candidate mRNA/lncRNA
  #' @param		cnv_M:		copy number profiles after filter
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		ceRNA_mir_list:		The interactions of mRNA/lncRNA-miRNA 
  #' @param		mir_ceRNA_list:		The interactions of miRNA-mRNA/lncRNA
  #' @param		cthr1:		Correlation coefficient threshold,mir与ceRNA1
  #' @param		cthr2:		Correlation coefficient threshold,mir与ceRNA2
  #' @param		cthr3:		Correlation coefficient threshold,ceRNA2与ceRNA1
  #' @param		pthr:		Correlation coefficient significance threshold
  dir.create(geneid,mode = "775")###create a folder for each gene, save the middle results
  temp_path<-getwd();
  temp_path<-paste(temp_path,"/",geneid,sep="")
  setwd(temp_path)##
  ceRNA_triplets<-ceRNA_mir_ceRNA(geneid,ceRNA_mir_list,mir_ceRNA_list)##ceRNA,miRNA,ceRNA
  save(ceRNA_triplets,file="ceRNA_triplets.Rdata")
  flag_v<-cnv_M[geneid,]
  exprs<-case_normal_profile(flag_v,gene_exp,mir_M) ##
  cnv_flag<-names(exprs)
  ##Significantly ceRNA triplets for each state
  sig_ceRNAs<-lapply(cnv_flag,cor_alls,exprs=exprs,ceRNA_triplets=ceRNA_triplets,ceRNA_mir_list=ceRNA_mir_list,mir_ceRNA_list=mir_ceRNA_list,cthr1=cthr1,cthr2=cthr2,cthr3=cthr3,pthr=pthr)
  ##
  names(sig_ceRNAs)<-cnv_flag
  save(sig_ceRNAs,file="sig_ceRNAs.Rdata")
  setwd("../")
  return(sig_ceRNAs)
  
  
}
###################################################################
##identify ceRNAs triplets
driver_lnc<-function(Dname,input,cthr1,cthr2,cthr3){
  #' @param		Dname:		disease name
  #' @param		input:		Folder path
  #' @param		cthr1:		Correlation coefficient threshold,mir与ceRNA1
  #' @param		cthr2:		Correlation coefficient threshold,mir与ceRNA2
  #' @param		cthr3:		Correlation coefficient threshold,ceRNA2与ceRNA1
  #' @param		p_thr:		Correlation coefficient significance threshold
  
  
  ##
  path_save<-paste(input,"/",Dname,sep="")
  setwd(path_save)
  #load("s1_result.RData")
  ## load("ce_mir.RData")
  ## cnv_M<-cnv_M[sig_genes,]
  
  #####
  save(cnv_M,file="cnv_M_for_driver.Rdata")
  save(cnv_M,gene_exp,mir_M,file = "filter_data.Rdata")
  genes<-rownames(cnv_M)
  result<-lapply(genes,ceRNAs_sig,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp,ceRNA_mir_list=ceRNA_mir_list,mir_ceRNA_list=mir_ceRNA_list,cthr1=cthr1,cthr2=cthr2,cthr3=cthr3,pthr=p_thr)
  names(result)<-genes
  return(result)
  
}

#########################################################################################################################

#extract ceRNA triplets(not null )
is.null<-function(a_list){
  
  if(length(unlist(a_list))==0){
    return(T)
  }else{
    return(F)
  }
  
}
extract_notnull<-function(all_sig_ceRNAs){
  #' @param		all_sig_ceRNAs:		genes affected by copy number and their triplets(result of driver_lnc)
  T_F<-unlist(lapply(all_sig_ceRNAs,is.null))
  result<-all_sig_ceRNAs[!T_F]
  return(result)
}

##################################################################################
##remove ceRNA triplets which ceRNA1 is the same as ceRNA2
remove_self1<-function(ceRNAs){
  if(length(ceRNAs)==0){
    return(ceRNAs)
  }else{
    if(length(ceRNAs)==3){
      return(NULL)
    }
    
    pos<-which(ceRNAs[,1]==ceRNAs[,3])
    if(length(pos)){
      ceRNAs<-ceRNAs[-pos,]
    }
    return(ceRNAs)
  }
}

remove_self2<-function(a_list){
  a_list<-lapply(a_list,remove_self1)
  return(a_list)
}
remove_self<-function(a_list){
#' @param		all_sig_ceRNAs:		genes affected by copy number and their triplets(result of driver_lnc and extract_notnull) 
  a_list<-lapply(a_list,remove_self2)
  return(a_list)
}

##############
#### rename ceRNA triplets 
rename_all_sig_ceRNAs<-function(all_sig_ceRNAs){
  #' @param		all_sig_ceRNAs:		genes affected by copy number and their triplets
  re<-lapply(all_sig_ceRNAs, function(i){
    temp<-lapply(i, function(j){
      if(is.null(j)){
        return(NULL)
      }else if(length(j)==3){
        ##	If the triple contains only three elements, convert to a matrix and concatenate the three elements (miRNA, mRNA or lncRNA) as the rowname
        j<-matrix(j,nrow = 1)
        rownames(j)<-paste(j[,1],j[,2],j[,3],sep = "_")
        return(j)
      }else{
        rownames(j)<-paste(j[,1],j[,2],j[,3],sep = "_")
        return(j)
      }
      
    })
    return(temp)
  })
  return(re)
}


###################
###	get ceRNA triplets and its ΔR
get_dysregulated_cmc<-function(all_sig_ceRNAs){
  ##all_sig_ceRNAs   genes affected by copy number and their triplets
  load("s1_result.RData")
  genes<-names(all_sig_ceRNAs)
  result<-lapply(genes,dysregulated_cmc,all_sig_ceRNAs,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp)
  names(result)<-genes
  return(result)
}


#############
##	for a gene affected by copy number,identify dysregulation extent of its triplets  
dysregulated_cmc<-function(geneid,all_sig_ceRNAs,cnv_M,mir_M,gene_exp){
  #' @param		geneid:		name(gene symbol) of this gene
  #' @param		all_sig_ceRNAs:		genes affected by copy number and their triplets
  #' @param		cnv_M:		copy number profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  flag_v<-cnv_M[geneid,]
  exprs<-case_normal_profile(flag_v,gene_exp,mir_M)
  cnv_flag<-names(exprs)
  ##	Calculate correlation coefficients in different states and compare 
  dysregulated<-lapply(cnv_flag, function(i){
    ceRNA_triplets<-all_sig_ceRNAs[[geneid]][[i]]
    if(is.null(ceRNA_triplets)){
      return(NULL)
    }else{
      sig_ceRNAs_c_p<-lapply(cnv_flag,get_cors,exprs=exprs,ceRNA_triplets=ceRNA_triplets)
      names(sig_ceRNAs_c_p)<-cnv_flag
      R<-abs(sig_ceRNAs_c_p[[1]][,3]-sig_ceRNAs_c_p[[2]][,3])## ΔR
      if(length(R)>1){
        #pos<-which(R>0.2)
        r<-R
        cmc<-ceRNA_triplets
        if(length(cmc)==0){
          result<-NULL
        }else if(length(cmc)==3){
          cmc<-matrix(cmc,nrow = 1)
          rownames(cmc)<-paste(cmc[,1],cmc[,2],cmc[,3],sep = "_")
          names(r)<-paste(cmc[,1],cmc[,2],cmc[,3],sep = "_")
          result<-list(cmc,r)
          names(result)<-c("dysregulated_cmc","R")
        }else{
          names(r)<-rownames(cmc)
          result<-list(cmc,r)
          names(result)<-c("dysregulated_cmc","R")
        }
        
      }else if(length(R)==1){
        #if(R>0.5){
          cmc<-ceRNA_triplets
          r<-R
          rownames(cmc)<-paste(cmc,collapse = "_")
          names(r)<-paste(cmc,collapse = "_")
          result<-list(cmc,r)
          names(result)<-c("dysregulated_cmc","R")
        #}else{
          #result<-NULL
        #}
      }
      
      return(result)
    }
  })
  
  # print(geneid)
  names(dysregulated)<-cnv_flag
  return(dysregulated)
}



############
#####get_cors calls get_cor,return Correlation coefficient between ceRNA1 and ceRNA2
get_cors<-function(a_flag,exprs,ceRNA_triplets){
  #' @param		a_flag:		cnv status tag
  #' @param		exprs:		expression of miRNA,lnc/mRNA
  #' @param		ceRNA_triplets:		ceRNA triplets in the current state
  a_list<-exprs[[a_flag]]
  gene_exp<-a_list[[1]]
  
  mir_M<-a_list[[2]]
  result<-get_cor(ceRNA_triplets,mir_M,gene_exp)
  return(result)
}

#######
get_cor<-function(ceRNA_triplets,mir_M,gene_exp){
  #' @param		ceRNA_triplets:		ceRNA triplets in the current state
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  if(length(ceRNA_triplets)==3){
    ceRNA1_ceRNA2<-c(ceRNA_triplets[1],ceRNA_triplets[3])
    ceRNA1_ceRNA2<-matrix(ceRNA1_ceRNA2,nrow = 1)
  }else{
    ceRNA1_ceRNA2<-ceRNA_triplets[,c(1,3)]
  }
  
  colnames(ceRNA1_ceRNA2)<-c("ceRNA1","ceRNA2")
  
  ceRNA1_ceRNA2_c_p<-only_corrs(ceRNA1_ceRNA2,expr1=gene_exp,expr2=gene_exp)
  return(ceRNA1_ceRNA2_c_p)
}

####### caculate Correlation coefficient between ceRNA1 and ceRNA2,only_corrs calls only_cors
only_cors<-function(a_pair,expr1,expr2){
  #' @param		a_pair:		a vector with 2 RNA names 
  #' @param		expr1和expr2:		Representing the corresponding expression profile
  ##print(a_pair)
  temp<-cor.test(expr1[a_pair[1],],expr2[a_pair[2],])
  pcor<-temp$estimate
  p<-temp$ p.value
  result<-c(pcor,p)
  return(result)
  
}
######
only_corrs<-function(all_pairs,expr1,expr2){
  #' @param		all_pairs:		all RNA names pairs
  #' @param		expr1和expr2:		Representing the corresponding expression profile
  c_p<-t(apply(all_pairs,1,only_cors,expr1=expr1,expr2=expr2))
  colnames(c_p)<-c("cor","pvalue")
  c_p<-cbind(all_pairs,as.data.frame(c_p))
  return(c_p)
}


########
#####扰乱
random_permute_filter<-function(all_dysregulated_cmc){
  #' @param		all_dysregulated_cmc:		genes affected by copy number and their triplets
  ##load("s1_result.RData")
  genes<-names(all_dysregulated_cmc)
  all_sig_dysregulated_cmc<-lapply(genes,random_permute_tests,all_dysregulated_cmc,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp)
  names(all_sig_dysregulated_cmc)<-genes
  return(all_sig_dysregulated_cmc)
}
###########
###random_permute_tests calls random_permute_test
random_permute_test<-function(geneid,all_dysregulated_cmc,cnv_M,mir_M,gene_exp){
  #' @param		all_dysregulated_cmc:		genes affected by copy number and their triplets
  #' @param		geneid:		name(gene symbol) of this gene
  #' @param		cnv_M:		copy number profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  flag_v<-cnv_M[geneid,]
  s_names<-sample(names(flag_v),length(flag_v))###randomize the labels
  names(flag_v)<-s_names
  exprs<-case_normal_profile(flag_v,gene_exp,mir_M)
  cnv_flag<-names(exprs)
  ##return R after randomize the labels
  dysregulated<-lapply(cnv_flag, function(i){
    ceRNA_triplets<-all_dysregulated_cmc[[geneid]][[i]][[1]]
    # print(1)
    if(is.null(ceRNA_triplets)){
      return(NULL)
    }else{
      sig_ceRNAs_c_p<-lapply(cnv_flag,get_cors,exprs=exprs,ceRNA_triplets=ceRNA_triplets)
      names(sig_ceRNAs_c_p)<-cnv_flag
      R<-abs(sig_ceRNAs_c_p[[1]][,3]-sig_ceRNAs_c_p[[2]][,3])
      if(length(R)>1){
        names(R)<-rownames(ceRNA_triplets)
        result<-R
      }else if(length(R)==1){
        names(R)<-paste(ceRNA_triplets,collapse = "_")
        result<-R
      }
      
      return(result)
    }
    # print(i)
  })
  
  names(dysregulated)<-cnv_flag
  return(dysregulated)
}
#################
###
random_permute_tests<-function(geneid,all_dysregulated_cmc,cnv_M,mir_M,gene_exp){
  #' @param		all_dysregulated_cmc:		genes affected by copy number and their triplets
  #' @param		geneid:		name(gene symbol) of this gene
  #' @param		cnv_M:		copy number profiles after filter
  #' @param		mir_M:		miRNA expression profiles after filter
  #' @param		gene_exp:		mRNA/lncRNA expression profiles after filter
  re<-lapply(1:1000, function(i){
    permute<-random_permute_test(geneid,all_dysregulated_cmc,cnv_M,mir_M,gene_exp)
    return(permute)
  })
  len<-lengths(re)[1]
  ##The randomized R is integrated into a matrix,row is triplets,col is R for each loop
  step_1<-lapply(1:len,function(i){## Two states
    step_2<-lapply(1:length(re), function(j){## 1000 loop
      step_3<-re[[j]][[i]]
    })
    if(is.null(step_2[[1]])){
      result<-NULL
    }else{
      result<-do.call("cbind",step_2)## matrix
    }
    return(result)
  })
  # print(dim(step_1[[1]]))
  
  # print(dim(step_1[[2]]))
  # print(rownames(step_1[[2]]))
  # print(geneid)
  names(step_1)<-names(re[[1]])
  ##Calculate p-values for randomly disturbed R matrix
  sig_dysregulated_cmc<-sig_p(geneid = geneid,step_1,all_dysregulated_cmc = all_dysregulated_cmc)
  # print(geneid)
  return(sig_dysregulated_cmc)
}

sig_p<-function(geneid,a_list,all_dysregulated_cmc){
  ##Calculate p-values for randomly disturbed R matrix
  #' @param		geneid:		name(gene symbol) of this gene
  #' @param		a_list:		R matrix after randomize 1000 times
  #' @param		all_dysregulated_cmc:		genes affected by copy number and their triplets
  cnv_flag<-names(a_list)
  step_1<-lapply(cnv_flag, function(i){
    real_R<-all_dysregulated_cmc[[geneid]][[i]][[2]]##  R under real conditions
    # print(real_R)
    if(is.null(a_list[[i]])){
      result<-NULL
    }else if(length(a_list[[i]])==1000){##  if only one triplets
      temp<-matrix(a_list[[i]],nrow = 1)
      re<-apply(temp,2,function(j){## Compare whether R of each loop is greater than R which under real conditions
        is_greater<-(j>real_R)
        return(is_greater)
      })
      p<-sum(re)/1000## Calculate the p value
      # print(p)
      if(p<0.05){## return significant dysregulated triplets
        result<-all_dysregulated_cmc[[geneid]][[i]][[1]]
        result<-matrix(unlist(result),nrow = 1)
        rownames(result)<-paste(result[,1],result[,2],result[,3],sep = "_")
      }else{
        result<-NULL
      }
    }else{
      re<-apply(a_list[[i]],2,function(j){
        is_greater<-(j>real_R)
        return(is_greater)
      })
      # print(re)
      p<-rowSums(re)/1000
      pos<-which(p<0.05)
      result<-all_dysregulated_cmc[[geneid]][[i]][[1]][pos,]
      if(length(result)==0){
        result<-NULL
      }else if(length(result)==3){
        result<-matrix(result,nrow = 1)
        rownames(result)<-paste(result[,1],result[,2],result[,3],sep = "_")
      }
    }
    
    return(result)
  })
  names(step_1)<-cnv_flag
  return(step_1)
  
}

########R under real conditions of triplets
get_real_R<-function(all_dysregulated_cmc,all_sig_dysregulated_cmc){
  #' @param		all_dysregulated_cmc:		result of get_dysregulated_cmc
  #' @param		all_sig_dysregulated_cmc:		result of random_permute_filter
  genes<-names(all_dysregulated_cmc)
  result<-lapply(genes,get_r,all_dysregulated_cmc,all_sig_dysregulated_cmc)
  names(result)<-genes
  return(result)
}

###
get_r<-function(geneid,all_dysregulated_cmc,all_sig_dysregulated_cmc){
  #' @param		geneid:		name(gene symbol) of this gene
  #' @param		all_dysregulated_cmc:		result of get_dysregulated_cmc
  #' @param		all_sig_dysregulated_cmc:		result of random_permute_filter
  cnv_flag<-names(all_dysregulated_cmc[[geneid]])
  step_1<-lapply(cnv_flag, function(i){
    sig_cmcs<-all_sig_dysregulated_cmc[[geneid]][[i]]
    if(length(sig_cmcs)==0){
      result<-NULL
    }else if(length(sig_cmcs)==3){
      tag<-paste(sig_cmcs,collapse = "_")
      r<-all_dysregulated_cmc[[geneid]][[i]][[2]][tag]
      result<-r
    }else{
      tag<-rownames(sig_cmcs)
      r<-all_dysregulated_cmc[[geneid]][[i]][[2]][tag]
      result<-r
    }
    return(result)
  })
  names(step_1)<-cnv_flag
  return(step_1)
}
########

##########get triplets with information of up or down 
get_up_down_cmc<-function(all_sig_dysregulated_cmc){
#' @param		all_sig_dysregulated_cmc:		result of random_permute_filter
  load("s1_result.RData")
  genes<-names(all_sig_dysregulated_cmc)
  result<-lapply(genes,up_down_cmc,all_sig_dysregulated_cmc,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp)
  names(result)<-genes
  return(result)
}

####
# decide the triplet is up or down
up_down_cmc<-function(geneid,all_sig_dysregulated_cmc,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp){
  flag_v<-cnv_M[geneid,]
  exprs<-case_normal_profile(flag_v,gene_exp,mir_M)
  cnv_flag<-names(exprs)
  ##
  u_d_cmc<-lapply(cnv_flag, function(i){
    ceRNA_triplets<-all_sig_dysregulated_cmc[[geneid]][[i]]
    if(is.null(ceRNA_triplets)){
      return(NULL)
    }else{
      sig_ceRNAs_c_p<-lapply(cnv_flag,get_cors,exprs=exprs,ceRNA_triplets=ceRNA_triplets)
      names(sig_ceRNAs_c_p)<-cnv_flag
      # print(paste(cnv_flag,collapse = ""))
      #cor.v(ceRNA1, ceRNA2) – cor.n(ceRNA1, ceRNA2)
      if(paste(cnv_flag,collapse = "")=="-20"){
        difference<-sig_ceRNAs_c_p[[1]][,3]-sig_ceRNAs_c_p[[2]][,3]
      }else{
        difference<-sig_ceRNAs_c_p[[2]][,3]-sig_ceRNAs_c_p[[1]][,3]
      }
      # difference<-sig_ceRNAs_c_p[[1]][,3]-sig_ceRNAs_c_p[[2]][,3]
      if(length(difference)>1){
        x<-c()
        x[which(difference>=0)]<-"up"
        x[which(difference<0)]<-"down"
        result<-cbind(x,ceRNA_triplets)
        
        
      }else if(length(difference)==1){
        if(difference>=0){
          x<-"up"
        }else{
          x<-"down"
        }
        ceRNA_triplets<-matrix(ceRNA_triplets,nrow=1)
        rownames(ceRNA_triplets)<-paste(ceRNA_triplets,collapse = "_")
        result<-cbind(x,ceRNA_triplets)
        
      }
      
      return(result)
    }
  })
  names(u_d_cmc)<-cnv_flag
  return(u_d_cmc)
}
##########

##########



# Merging triplets of two states
get_union_cmc<-function(all_sig_dysregulated_cmc_up_down){
  #' @param		all_sig_dysregulated_cmc_up_down:		result of get_up_down_cmc
  # genes<-names(all_sig_dysregulated_cmc_up_down)
  result<-lapply(all_sig_dysregulated_cmc_up_down,function(i){
    re<-do.call("rbind",i)
    re<-unique(re)
    return(re)
  })
  # names(result)<-genes
  return(result)
}

