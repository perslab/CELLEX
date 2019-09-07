import numpy as np
import pandas as pd
import plotnine as p9

def n_es_genes(df: pd.DataFrame, 
                annotation: pd.Series,
                figsize: tuple=None) -> p9.ggplot:
    """Plot distribution of number of ES genes per group
    
    Computes the number of ES genes per column, e.g. cell(-type) 
    and plots the distribution for the groups specified
    by the annotation.
    
    Parameters
    ----------
    df : DataFrame
        Dataframe containing positive ES weights, ideally use only ESmu.
    annotation : Series
        Annotation to group dataframe cell(-types) by in the violin plots.
    figsize : (float, float), optional (default: None)
        Specify width and height of plot.
    
    Returns
    -------
    p : ggplot
        A plotnine ggplot

    """
    
    ### Count number of non-zero values, i.e. ESw > 0
    df = df.astype(bool).sum(axis=0)
    
    ### Map column labels to annotation
    if type(annotation) is pd.DataFrame:
        annotation = annotation.iloc[:,0]
    
    # remove duplicates
    annotation = annotation.loc[~annotation.index.duplicated(keep='first')]

    df.index = df.index.map(annotation, na_action="ignore").values.astype(str)
    
    # Constants, height and width of plot.
    if figsize is None:
        W = min((df.index.nunique(), 10))
        H = 6.4 # plotnine default height
    else:
        W, H = figsize

    ### Convert to tidy / long format
    # Org:
    #       ABC  ACBG  ACMB
    # POMC  0.0   0.5   0.9
    # AGRP  0.2   0.0   0.0
    # LEPR  0.1   0.1   0.4
    
    # Tidy:
    #   gene_name annotation    es_weight
    # 1 POMC      ABC           0.0
    # 2 AGRP      ABC           0.6
    # 3 LEPR      ABC           1.0      
    df_tidy = df.copy()
    df_tidy.index.name = None
    df_tidy = pd.melt(df_tidy.reset_index(), id_vars="index", var_name="annotation", value_name="count")
    
    ### Compute the mean count of ES genes
    mean_count = df_tidy["count"].mean(axis=0)
    
    ### Plot
    p = (
        ### data
        p9.ggplot(data=df_tidy, 
                  mapping=p9.aes(x="index", y="count", fill="index", label="index"), 
                 )

        ### theming
        + p9.theme_classic()
        + p9.theme(
            figure_size = (W,H),
            axis_text_x = p9.element_text(rotation=75)
        )

        + p9.labs(
            x="", # e.g. "Cell-type"
            y="Number of ES genes", # e.g. "ES weight"
        )
        
        ### viz
        + p9.geom_violin(scale="width", show_legend=False)
        + p9.geom_jitter(width=0.1, height=0, show_legend=False)
        + p9.geom_hline(yintercept=mean_count, color="blue", linetype="dashed", show_legend=False)
    )
    
    return p
