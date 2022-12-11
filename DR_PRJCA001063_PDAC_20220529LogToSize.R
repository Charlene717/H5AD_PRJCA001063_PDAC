##### Current path and new folder setting* #####
ProjectName = "TOP2A"
Sampletype = "PDAC"
#ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){
  dir.create(Save.Path)
}

IGene = "TOP2A"



# cds_sub_AcinaDucT_NewK_ReCluster_TOP2A
# cds_sub_DucT2_TOP2A_SmallCTR_TOP2A

cds_sub_DucT2MDC <- choose_cells(cds)

Main = IGene


plot_cells(cds_sub_DucT2MDC , genes = Main, show_trajectory_graph = F,label_cell_groups = F
           ,cell_size =1,norm_method = c("size_only") )+
  scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  ggtitle(Main)+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, 
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 0.8),
        axis.line.y = element_line(colour = "black", size = 0.8)) -> P

P


pdf(
  file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_cds_sub_DucT2MDC_",IGene,".pdf"),
  width = 7, height = 7
)

print(P)

graphics.off()

