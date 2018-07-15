function () 
{
	library(gplots)
	library(xlsx)

	options(stringsAsFactors=FALSE)
	#par(ask=T)

	my_data=read.csv("data/paired_tpm_20",sep="\t",header=T,row.names=1)
	my_data_scaled=scale(my_data)
	xx=unlist(strsplit(rownames(my_data),"\\."))
	xx=xx[1+2*(0:((length(xx)-1)/2))]
	my_data_scaled=cbind(my_data_scaled,xx)
	
	#	---------------------------------------------------
	#
	#	heatmaps depicting the clustering of the samples
	#
	#	---------------------------------------------------
	heatmap_sample_clustering=function(){
		heatmap.2(cor(log(my_data+1)))	
		heatmap.2(cor(my_data_scaled))	
		B=c()
		for(x in 1:dim(my_data)[2]){
			D=subset(my_data,my_data[,x]>1)
			cl=MAP[rownames(D)]
			B=c(B,length(unique(cl)))
		}
		barplot(B,names=colnames(my_data),las=2,ylim=c(0,25000))
	}

	#
	#	orthologous screen depicting the Arabidopsis-Potato interaction
	#
	compare_to_jens_data=function(){
		data=read.xlsx("data/orthologous_screen.xlsx",1)	
		xx=subset(data,data[,6]=="2")[,2]
		my_goi=subset(my_data_scaled,my_data_scaled[,dim(my_data_scaled)[2]] %in% xx)
		my_goi=my_goi[,2:(dim(my_goi)[2]-1)]
		my_goi2=matrix(as.double(my_goi),nrow=dim(my_goi)[1])
		print((my_goi2))
		pheatmap(my_goi2)
	}		

	#	---------------------------------------------------
	#
	# heatmap sample clustering
	#
	#	---------------------------------------------------

	heatmap_sample_clustering()
	compare_to_jens_data()
}
