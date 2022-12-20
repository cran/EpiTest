#' Display variance component graphs
#'
#' @param ETest a list as obtained from the EpiTest.fit function
#' @param Title (optional) a character string to be used as a title for the graph
#' @param Family_names (optional) a character vector with family names
#' @param Colors (optional) a character vector with three colors for the
#' graphical display
#'
#' @import ggplot2
#' @import stringr
#' @import dplyr
#'
#' @return a ggplot graph
EpiTest.plot.var <- function(ETest,Title = NULL,Family_names = NULL,Colors = NULL){

  ## Vector of family names
  Families <- names(ETest) %>% set_names

  ## dataframe with variances and associated pval according to family and component
  Var_df <- Families %>%
    map(~ETest[[.x]]$Sigma2) %>%
    imap_dfr(~tibble(Family = .y,Component = names(.x) %>% stringr::str_split("_") %>% map_chr(last),Variance = .x)) %>%
    left_join(ETest %>% imap_dfr(~c(Pval=.x$Test_random$pval[2],Family=.y,Component="SxS")),
              by = c("Family", "Component")) %>%
    mutate(Test=ifelse(.data$Component%in%c("S","E"),"",
                       ifelse(.data$Pval<0.001,"***",
                              ifelse(.data$Pval<0.01,"**",
                                     ifelse(.data$Pval<0.05,"*",""))))) %>%
    # mutate(Variance=ifelse(.data$Component=="S",.data$Variance*0.25,
    #                        ifelse(.data$Component=="SxS",.data$Variance*0.0625,.data$Variance))) %>%
    mutate(Test=as.factor(.data$Test)) %>%
    mutate(Component=factor(.data$Component,levels=c("E","SxS","S")))

  ## Change family names
  if(!is.null(Family_names)){
    Var_df <- Var_df %>%
      left_join(bind_cols(Family=Families,Family_new=Family_names),by = "Family") %>%
      mutate(Family=.data$Family_new)
  }

  ## Colors
  if(is.null(Colors)){
    Colors <- c(E = "grey",S = "steelblue2",SxS = "orange")
  }else{
    Colors <- Colors %>% set_names(nm = c("E","S","SxS"))
  }

  ## Plot
  EpiVarPlot <- ggplot2::ggplot(Var_df,aes(x=.data$Family,y=.data$Variance,fill=.data$Component)) +
    theme_bw() +
    labs(title = Title) +
    xlab("Population")+
    scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    scale_fill_manual(values = Colors,labels=c(E=expression(sigma[E]^2),
                                               S=expression(sigma[S]^2),
                                               SxS=expression(sigma[SxS]^2))) +
    # scale_fill_manual(values = Colors,labels=c(E=expression(sigma[E]^2),
    #                                            S=expression(pi[i]*(1-pi[i])*sigma[S]^2),
    #                                            SxS=expression((pi[i]*(1-pi[i]))^2*sigma[SxS]^2))) +
    geom_bar(position="stack",stat="identity") +
    geom_text(aes(label=.data$Test),position=position_stack(vjust=0.5))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text.align = 0,
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.margin = margin(r=20))

  return(EpiVarPlot)
}


#' Display directional epistasis plot
#'
#' @param ETest a list as obtained from the EpiTest.fit function
#' @param Title (optional) a character string to be used as a title for the graph
#' @param Family_names (optional) a character vector with family names
#' @param Alpha type I error for the test
#'
#' @return A ggplot graph
EpiTest.plot.mean <- function(ETest,Title = NULL,Family_names = NULL,Alpha = NULL){

  ## p-value
  Pval_delta <- ETest %>% map_dbl(~.x$Test_fixed[3,2])

  ## Dataframe with population parameters
  Param_df <- ETest %>%
    imap_dfr(~c(Family=.y,.x$Beta)) %>%
    mutate_at(c("(Intercept)","beta","delta"),as.numeric) %>%
    mutate(Pval_delta = Pval_delta) %>%
    mutate(Size = ifelse(.data$Pval_delta>Alpha,0.25,0.75))

  ## dataframe with parents parameters
  Parent_df <- tibble(x = c(0,rep(1,nrow(Param_df))),
                      y = c(Param_df[,"(Intercept)"] %>% colMeans,
                            Param_df[,c("(Intercept)","beta","delta")] %>% rowSums),
                      Names = c("",Param_df$Family))

  ## Change family names
  if(is.null(Family_names)==F){
    Parent_df$Names[-1] <- Family_names %>% as.character
  }

  ## Plot
  expand_xleft <- (nchar(Parent_df$Names) %>% max)*0.017 + 0.03
  EpiDirPlot <- ggplot(Parent_df,aes(.data$x,.data$y,label=.data$Names)) +
    labs(title = Title)+
    xlab("Alternative parent proportions") +
    ylab("Expected response") +
    scale_x_continuous(expand = expansion(mult = c(0.03,expand_xleft))) +
    apply(Param_df,1,function(z)stat_function(
      fun = function(x, `(Intercept)`, beta, delta) {`(Intercept)` + x*beta + x^2*delta},
      args = list(`(Intercept)` = z["(Intercept)"] %>% as.numeric,
                  beta = z["beta"] %>% as.numeric,
                  delta = z["delta"] %>% as.numeric),
      aes(color = z["Pval_delta"] %>% as.numeric %>% log10 %>% abs),
      size = z["Size"] %>% as.numeric)) +
    scale_color_gradient2(low = "gold",mid = "red",high = "red4",midpoint = 12,na.value = "grey",
                          limits=c(-log10(Alpha),max(c(-log10(min(Param_df$Pval_delta)),20))+1),
                          breaks = seq(0,max(c(-log10(min(Param_df$Pval_delta)),20))+1,5),
                          name=expression("Test"~delta~(-log[10](italic(p))))) +
    ggrepel::geom_text_repel(xlim = c(1,Inf),direction = "y",size=2.5,max.overlaps = 20,segment.color="grey") +
    theme_bw()
  return(EpiDirPlot)
}


#' Display variance component or directional epistasis graphs
#'
#' @param ETest a list as obtained from the EpiTest.fit function
#' @param Title (optional) a character string to be used as a title for the graph
#' @param Family_names (optional) a character vector with family names
#' @param Colors (optional) a character vector with three colors for the
#' graphical display (for variance graphs only)
#' @param Type a string with "mean" to display the directional epistasis plot or
#' "var" to display the variance component graph
#' @param Alpha type I error for the test, with 5% default value
#' (for directional epistasis plot only)
#'
#' @return a ggplot graph
#' @export
#'
#' @examples
#' data(Pheno.list)
#' data(Ancestry.list)
#'
#' ## Full NAM analysis
#' Parent.list <- Parents$Genotype[-1]
#' names(Parent.list) <- Parents$Family[-1]
#' ETest.nam <- purrr::imap(Parent.list, ~ EpiTest.fit(Ancestry = Ancestry.list[[.y]],
#'                                                     Pheno = Pheno.list[[.y]],
#'                                                     ParentName = c(P0 = 'B73',P1 = .x),
#'                                                     Parents = Parents,
#'                                                     Trait = 'GDD_DTA',
#'                                                     Weight = TRUE))
#' ## Variance component plot
#' EpiTest.plot(ETest.nam,Title = 'Days to anthesis',Type = "var")
#'
#' ## Directional epistasis plot
#' EpiTest.plot(ETest.nam,Title = 'Days to anthesis',Type = "mean",Alpha = 5/100)
EpiTest.plot <- function(ETest,Title = NULL,Family_names = NULL,Colors = NULL,Type = "mean",Alpha = 5/100){
  if(Type=='var'){
    Res <- EpiTest.plot.var(ETest,Title = Title,Family_names = Family_names,Colors = Colors)
  }
  if(Type=='mean'){
    Res <- EpiTest.plot.mean(ETest,Title = Title,Family_names = Family_names,Alpha = Alpha)
  }
  return(Res)
}
