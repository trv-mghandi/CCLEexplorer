function(input, output, session) {
    
  #output$data = data.frame(feature = sample(colnames(F),4), F.type=sample(F.type, 4), r=sample(100,4)/100)
  
   output$plot <- renderPlotly({
      plot_ly(x = rnorm(10), y = rnorm(10),  type = "scatter")
    })

   observe({
     x <- f.types.id[match(input$query.feature.type, f.types.nam)]
     
     updateSelectInput(session, "query.feature.name",
                       label = "feature",
                       choices = colnames(F)[F.type==x]
     )
   })
   
    observe({
      x <- input$x.dataset
      
      updateSelectInput(session, "x.gene",
                        label = "Choose a feature:", 
                        choices = colnames(F)[F.type==x]#choices = g$Data[[x]]$name
      )
    })
    
    observe({
      y <- input$y.dataset
      
      updateSelectInput(session, "y.gene",
                        label = "Choose a feature:", 
                        choices = colnames(F)[F.type==y]#choices = g$Data[[y]]$name
                        
      )
    })
    
    observe({
      z <- input$z.dataset
      
      updateSelectInput(session, "z.gene",
                        label = "Choose a feature:", 
                        choices = colnames(F)[F.type==z]#choices = g$Data[[z]]$name
      )
    })
    
    observeEvent(input$check.all,{
      updateCheckboxGroupInput(session,"Filt.types", selected=f.types.id); 
      
    })
    observeEvent(input$uncheck.all,{
      updateCheckboxGroupInput(session,"Filt.types", selected=character(0)); 
      
    })
    
    observeEvent(input$refresh.plot,{
      output$plot <- renderPlotly({
        plot_ly(x = F[,isolate(input$x.gene)],#g$Data[[isolate(input$x.dataset)]]$dat[,isolate(input$x.gene)], 
                y = F[,isolate(input$y.gene)],#g$Data[[isolate(input$y.dataset)]]$dat[,isolate(input$y.gene)], 
                z = F[,isolate(input$z.gene)],#g$Data[[isolate(input$z.dataset)]]$dat[,isolate(input$z.gene)], 
                text = rownames(F),#g$Data[[isolate(input$x.dataset)]]$dat),
                #text = paste(rownames(g$Data[[isolate(input$x.dataset)]]$dat),'\n','x=',
                #             g$Data[[isolate(input$x.dataset)]]$dat[,isolate(input$x.gene)]),
                #hoverinfo = 'text',
                color = g$col, 
                colors = g$cols,  aspectmode='cube',
                type = ifelse(input$toggle.2D3D %%2==1,"scatter3d","scatter")) %>%         layout(scene = list(xaxis = list(title = paste0(isolate(input$x.gene))),
                                                                                                               yaxis = list(title = paste0(isolate(input$y.gene))),
                                                                                                               zaxis = list(title = paste0(isolate(input$z.gene)))))  %>%         layout(xaxis = list(title = paste0(isolate(input$x.gene))),
                                                                                                                                                                                                                       yaxis = list(title = paste0(isolate(input$y.gene))))
        
        
        #%>%
        #  add_markers() %>%
      })
      
    })

    observeEvent(input$toggle.all,{
      plotlyProxy("plot", session) %>%
        plotlyProxyInvoke("restyle", list(visible = ifelse(input$toggle.all%%2==1,'legendonly', TRUE)))
    })
    
    #input = list(); input$color.by='type_refined'
    observeEvent(input$color.by,{
      if(input$color.by %in% g$color_type_family)  {
        g$col<<-Annotations[,input$color.by] 
        coltab = subset(g$col_type_table, col_type_table$family==input$color.by)
        g$cols <<- coltab$col
        
        g$cols[which(is.na(g$cols))]<<-"#AAAAAA"
        g$col <<- factor(g$col, levels = c(coltab$type,'N/A'))
        g$col[which(is.na(g$col))]<<-'N/A'
        
        #g$cols= c(g$cols,"#AAAAAA" )
        
      }else{
        g$col<<-rep('all', nrow(g$Data[[1]]$dat))
        g$col<<- as.factor(g$col)
        g$cols<<- 'red'
      }
    })

    observe({
     #  idx = calcNetwork(geneName = 'ESR1',10,1,10,enets)
  # #    idx = calcNetwork(geneName = input$query.feature.name,as.numeric(as.character(input$query.nmax)),
  # #                    as.numeric(as.character(input$query.dmax)),as.numeric(as.character(input$query.emax)),enets)
     #  edges = edges_netformat[enets$edges_list[[idx]],]
    #   nodes = nodes_netformat[enets$nodes_list[[idx]],]

     idx = calcNetworkBigF(input$query.feature.name, nmax=as.numeric(as.character(input$query.nmax))+1, 
                           dmax=1,#as.numeric(as.character(input$query.dmax)), 
                           emax=as.numeric(as.character(input$query.emax)), 
                           f.type.ids.filt =input$Filt.types, enets)
#     f.type.ids.filt =g.cur.f.types.nam.selected, enets)
    
     edges = enets$edges_list[[idx]]
     nodes = enets$nodes_list[[idx]]
      
    #idx = calcNetwork(geneName = 'ESR1', 100000,1,100,enets)

     output$networkTopCodep2 <- renderVisNetwork({
        visNetwork(nodes, edges) %>%
          visEvents(select = "function(edges) {
                Shiny.onInputChange('selected_edges_ids', edges);
                ;}")  %>%  visEvents(doubleClick = "function(nodes){
                     Shiny.onInputChange('doubleClick', nodes.nodes[0]);
                      ;}")
       
       
       
       
      })

      #  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)%>%
      #  visGroups(groupname = "1", color = "maroon") %>% 

    })
    
    observe({
      #    input$gosel
      input$selected_edges_ids
      visNetworkProxy("networkTopCodep2") %>% visGetSelectedEdges()
      #print(paste(length(input$network_selectedEdges),input$network_selectedEdges))
    })
    
#     observe({
#       #    input$gosel
#       # input$doubleClick
#       #if(!is.null(input$doubleClick)){
#         ni = input$doubleClick; 
#         if(!is.null(ni)){
#           
#         if(ni!=0){
#           session$sendCustomMessage("doubleClick", 0)
#           ##x <- f.types.id[match(input$query.feature.type, f.types.nam)]
#           ti = F.type[ni]
#           fti = f.types.nam[match(ti,f.types.id)]
#           
# #          if(g$doubleClick == 0){ #to avoid two clicks bug 
# #            g$doubleClick = ni;
#             updateSelectInput(session, "query.feature.type",selected =fti);
#             updateSelectInput(session, "query.feature.name",selected =colnames(F)[ni]);
#  #         }else{
# #            updateSelectInput(session, "query.feature.name",selected =colnames(F)[ni]);
#  #           
# #            g$doubleClick = 0
#  #         }
#   
#       }}
#     })
    
    
    observe({
      #    input$gosel
      input$doubleClick
      if(!is.null(input$doubleClick)){
        ni = input$doubleClick; 
        
        if(ni!=0){
          session$sendCustomMessage("doubleClick", 0)
          x <- f.types.id[match(input$query.feature.type, f.types.nam)]
          ti = F.type[input$doubleClick]
          fti = f.types.nam[match(ti,f.types.id)]
          
          updateSelectInput(session, "query.feature.type",selected =fti);
          updateSelectInput(session, "query.feature.name",selected =colnames(F)[input$doubleClick]);
          
        }}
    })
    
    
    if(TRUE){  
      
    observe({
      input$networkTopCodep2_selectedEdges
      #      output$shiny_return <- renderPrint({paste(length(input$networkTopCodep2_selectedEdges),input$networkTopCodep2_selectedEdges)})
      g2 = '';
      if(length(input$networkTopCodep2_selectedEdges)==1){
        
       # qq=input$networkTopCodep2_selectedEdges
        
        #browser()
        elist =enets$edges_list[[enets$curIdx]]
        edgei = elist[input$networkTopCodep2_selectedEdges,];
        g1 = colnames(F)[edgei$from]
        g2 = colnames(F)[edgei$to]
        
        #g1 = genes[rsub.idx1[input$networkTopCodep2_selectedEdges]]
        #g2 = genes[rsub.idx2[input$networkTopCodep2_selectedEdges]]
        
        if(g2== isolate(input$query.feature.name)){
          hh = g1; g1=g2; g2=hh;
        }
        
      }  else{
        #qq=input$networkTopCodep2_selectedEdges   ##node selected 
        if(FALSE){
          g1 = isolate(input$query.feature.name)
          g2 = ""
          igi = which(genes==g1)
          jj = which(rsub.idx1[input$networkTopCodep2_selectedEdges]==igi); 
          if(length(jj)>0){
            g2 = genes[rsub.idx2[input$networkTopCodep2_selectedEdges[jj[1]]]]
          }else{
            jj = which(rsub.idx2[input$networkTopCodep2_selectedEdges]==igi); 
            if(length(jj)>0){
              g2 = genes[rsub.idx1[input$networkTopCodep2_selectedEdges[jj[1]]]]
            }            
          }
        }
        
      }
      if(g2!=''){
        
        #browser()
      
        ## replot 
        
        #updateSelectInput(session, "x.gene",selected =g1);
        #updateSelectInput(session, "y.gene",selected =g2);
        
        output$plot <- renderPlotly({
          plot_ly(x = F[,g1],#g$Data[[isolate(input$x.dataset)]]$dat[,isolate(input$x.gene)], 
                  y = F[,g2],#g$Data[[isolate(input$y.dataset)]]$dat[,isolate(input$y.gene)], 
                  z = F[,isolate(input$z.gene)],#g$Data[[isolate(input$z.dataset)]]$dat[,isolate(input$z.gene)], 
                  text = rownames(F),#g$Data[[isolate(input$x.dataset)]]$dat),
                  #text = paste(rownames(g$Data[[isolate(input$x.dataset)]]$dat),'\n','x=',
                  #             g$Data[[isolate(input$x.dataset)]]$dat[,isolate(input$x.gene)]),
                  #hoverinfo = 'text',
                  color = g$col, 
                  colors = g$cols,  aspectmode='cube',
                  type = ifelse(input$toggle.2D3D %%2==1,"scatter3d","scatter")) %>%         layout(scene = list(xaxis = list(title = paste0(g1)),
                                                                                                                 yaxis = list(title = paste0(g2)),
                                                                                                                 zaxis = list(title = paste0(isolate(input$z.gene)))))  %>%         layout(xaxis = list(title = paste0(g1)),
                                                                                                                                                                                           yaxis = list(title = paste0(g2)))

            })
        
        if(DEMO_mode){
           updateSelectInput(session, "x.dataset",selected =F.type[match(g1, colnames(F))]);  ## may make it slow
           updateSelectInput(session, "x.gene",selected =g1);
           updateSelectInput(session, "y.dataset",selected =F.type[match(g2, colnames(F))]);
           updateSelectInput(session, "y.gene",selected =g2);
        }
        
      }
      
      
       
    })
    
    
    }
     
#   output$shiny_return <- renderPrint({
#     #  paste(input$networkTopCodep_selected,' se')
# #    paste(input$query.feature.name,as.numeric(as.character(input$query.nmax))+0.1,as.numeric(as.character(input$query.dmax))+0.1,as.numeric(as.character(input$query.emax))+0.1,' se')
#     paste(input$query.feature.name,as.numeric(as.character(input$query.nmax)),as.numeric(as.character(input$query.dmax)),as.numeric(as.character(input$query.emax)),
#           match(as.numeric(as.character(input$query.nmax)),nmax.options),
#           match(as.numeric(as.character(input$query.dmax)),dmax.options),
#           match(as.numeric(as.character(input$query.emax)),emax.options))
#     #paste(' se')
#   })
    
}



















