import numpy as np
import pandas as pd
import plotnine as p9

def gene_profile(genes: list, 
                 weights: pd.DataFrame, 
                 stddev: pd.DataFrame=None,
                 y_axis_label: str=None,
                 highlight_n: int=None, 
                 highlight_anno: list=None, 
                 figsize: tuple=None,
                 ylim: tuple=None) -> p9.ggplot:
    """
    
    Parameters
    ----------
    weights            : DataFrame of ES weights
    genes          : a single str or list of genes to include in plot as facets
    highlight_n    : number of highest ESw to highlight
    highlight_anno : specific annotations to highlight
    figsize : (float, float), optional (default: None)
        Specify width and height of plot.
    
    Returns
    -------
        g    : ggplot
        
    Todo:
        * find a better way for sorting cell-types along x-axis
        * report if gene in genes is not found in df
        * report if duplicate genes
        * replace hacky x-axis labelling
    
    """
    
    ### Reduce dataframe to genes of interest
    genes = [str.upper(s) for s in genes]
    idx = np.char.upper(weights.index.values.astype(str))
    mask = np.isin(idx, genes)
    df_tidy = weights[mask]
    n_genes = len(df_tidy)

    assert (n_genes >= 1), "No matching genes found in dataframe."

    stddev_tidy = None
    if stddev is not None:
        idx = np.char.upper(stddev.index.values.astype(str))
        mask = np.isin(idx, genes)
        stddev_tidy = stddev[mask]
        n_genes = len(df_tidy)
        assert (n_genes >= 1), "No matching genes found in stddev dataframe."

    # Constants, height and width of plot.
    if figsize is None:
        H = 5*n_genes
        W = 15
    else:
        W, H = figsize

    if ylim is None:
        ylim = (-1,1)
    
    if y_axis_label is None:
        y_axis_label = "Expression Specificity"
    
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
    
    if stddev_tidy is not None:
        stddev_tidy.index.name = None
        stddev_tidy = pd.melt(stddev_tidy.reset_index(), id_vars="index", var_name="annotation", value_name="stddev")
        df_tidy = df_tidy.merge(stddev_tidy, on=["index", "annotation"])


    ### Sort values by gene_name and es_weight and add order
    # Sorted:
    #   gene_name annotation   es_weight   x_order
    # 1 AGRP      MOL2         0.0         1
    # 2 AGRP      ACNT1        0.1         2
    # 3 AGRP      MOL1         0.2         3
    
    df_tidy = df_tidy.sort_values(by=["index", "weight"])
    df_tidy["order"] = np.arange(len(df_tidy)) + 1
    
    ### Generate highlight
    # Default: highlight top 5
    if ((highlight_n is None) and (highlight_anno is None)):
        highlight_n = 5

    # highlight list of 
    if (highlight_anno is not None):
        df_tidy["highlight"] = df_tidy["annotation"].isin(highlight_anno)
    elif (highlight_n is not None):
        df_tidy["highlight"] = df_tidy.groupby("index")["order"].rank("first", ascending=False) <= highlight_n
    else:
        df_tidy["highlight"] = np.array([False] * len(df_tidy))
    
    df_highlight = df_tidy[df_tidy["highlight"]]
    
    ### Plot
    # linear function to compute x_axis text-size.
    # Mainly depends on number of genes in df per faceet, i.e. len(df_tidy) / len(genes).
    SIZE_TEXT_X_AXIS = 10.161 - 0.023 * (len(df_tidy) / len(genes))
    
    # Limits of the order for each index gene / facet, e.g. [0, 266, 531]
    # These limits are necessary to only plot the labels
    order_lims = [0, *(df_tidy.groupby("index")["order"].max().values)]
    
    def find_nearest(array,value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
        
    def getbreaks(lims):
        # function defined for use in debugging
        l = find_nearest(order_lims, lims[0])
        r = find_nearest(order_lims, lims[1])
        breaks = np.arange(l, r)
        return breaks

    def getlbls(idx):
        # function defined for use in debugging
        idx = idx
        lbls = df_tidy["annotation"].iloc[idx].values
        return lbls
    
    p = (
        ### data
        p9.ggplot(data=df_tidy, mapping=p9.aes(x="order", y="weight", label="annotation"))

        ### theming
        + p9.theme_classic()
        + p9.theme(
            figure_size = (W,H),
            axis_ticks_major_x = p9.element_blank(),
            axis_text_x = p9.element_text(rotation=75, hjust=0, size=SIZE_TEXT_X_AXIS), # 
            axis_text_y = p9.element_text(size=W),
            panel_spacing = 1,
            strip_background = p9.element_blank()
        )

        + p9.ylim(ylim[0],ylim[1])

        + p9.labs(
            x="", # e.g. "Cell-type"
            y=y_axis_label, # e.g. "ES weight"
        )

        ### viz
        # all
        + p9.geom_segment(mapping=p9.aes(x="order", xend="order", y=0, yend="weight"),
                       color="grey",
                       alpha=0.3,
                       show_legend=False
        )

        + p9.geom_point(mapping=p9.aes(size=2),
                     color="grey",
                    show_legend=False
        )

        # highlight
        + p9.geom_point(data=df_highlight, mapping=p9.aes(size=2), 
                     color="dodgerblue",
                    show_legend=False
        )

        + p9.geom_segment(data=df_highlight, mapping=p9.aes(x="order", xend="order", y=0, yend="weight"),
                       color="dodgerblue",
                       alpha=0.3,
                       show_legend=False
        )

        + p9.facet_wrap("index",
                     scales="free",
                     nrow=n_genes
                    )
        
        + p9.scale_x_continuous(
            # order_scale is continuous across all annotations
            # so the scale will look weird for each facet, e.g.
            # facet 1 may have order 1-7, and facet 2 has order 8-14.
            # therefore we must use a labeller function to get the 
            # correct labels for each interval of order.
            breaks = lambda lims: getbreaks(lims),
            labels = lambda idx: getlbls(idx)
        )
    )
    
    if stddev_tidy is not None:
        p = p + p9.geom_errorbar(mapping=p9.aes(ymin="weight-stddev", ymax="weight+stddev"), 
                                    color="grey", width=0.1)\
                + p9.geom_errorbar(data=df_highlight, mapping=p9.aes(ymin="weight-stddev", ymax="weight+stddev"),
                                color="dodgerblue", width=0.1)

    # add labels last for them to be on top
    p = p + p9.geom_label(data=df_highlight,
                    color = "dodgerblue",
                    adjust_text = {'expand_points': (2,2)}
        )

    return p
    