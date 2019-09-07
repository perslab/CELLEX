import numpy as np
import pandas as pd
import plotnine as p9

def heatmap(esw: pd.DataFrame, 
            genes: list=None, 
            annotations: list=None,
            figsize: tuple=None) -> p9.ggplot:
    """
    
    Args:
    esw             : DataFrame of ES weights
    genes          : a list of genes to include in the heatmap
    annotations    : a list of annotations to include in the heatmap
    figsize : (float, float), optional (default: None)
        Specify width and height of plot.
    Returns:
        g    : ggplot

    """

    df_tidy = esw

    ### Reduce dataframe to genes and annotations of interest
    if genes is not None:
        genes = [str.upper(s) for s in genes]
        idx = np.char.upper(df_tidy.index.values.astype(str))
        mask = np.isin(idx, genes)
        df_tidy = esw[mask]
    
    if annotations is not None:
        annotations = [str.upper(s) for s in annotations]
        cols = np.char.upper(df_tidy.columns.values.astype(str))
        mask = np.isin(cols, annotations)
        df_tidy = df_tidy.iloc[:,mask]
    
    # Constants, height and width of plot.
    if figsize is None:
        W = min((df_tidy.shape[0], df_tidy.shape[1], 20))
        H = min((df_tidy.shape[0], df_tidy.shape[1], 20))
    else:
        W, H = figsize

    ### Convert to tidy / long format if necessary
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

    df_tidy.index.name = None # ensure that index name is none, so "index" is used for id_vars
    df_tidy = pd.melt(df_tidy.reset_index(), id_vars="index", var_name="annotation", value_name="weight")
       
    ### Plot
    p = (
        ### data
        p9.ggplot(data=df_tidy, mapping=p9.aes(x="index", y="annotation", fill="weight", label="annotation"))

        ### theming
        + p9.theme_classic()
        + p9.theme(
            figure_size = (W,H),
            axis_text_x = p9.element_text(rotation=75),
        )

        + p9.labs(
            x="", # e.g. "Cell-type"
            y="", # e.g. "ES weight"
        )
        
        ### viz
        + p9.geom_tile()
        # + p9.scale_fill_gradientn(colors=['#9ebcda','#8c6bb1','#88419d','#6e016b']) # light blue to purple
        + p9.scale_fill_gradientn(colors=['#ffffff','#1E90FF'], limits=[0,1]) # white to dodgerblue
    )
    
    return p
