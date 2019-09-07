import os
import plotnine as p9
from datetime import datetime

def save_as_pdf(plot: p9.ggplot, filename: str=None, path: str=None, dpi: int=None, verbose: bool=False) -> None:
    """Save a plotnine ggplot as pdf

    Parameters
    ----------
    plot : p9.ggplot
        The plot to save

    filename : str, optional (default: None)
        Filename to write to. If None, a name is generated.

    path : str, optional (default: None)
        Path to save to. If None, saves to "out".

    dpi : int, optional (default: None)
        DPI of saved plot. If None, set to 300.

    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    None

    """
    
    if path is None:
        path = "out"
    
    if filename is None:
        dateTimeObj = datetime.now()
        filename = "{}/cellex_plot_{}.pdf".format(path, dateTimeObj.strftime("%y%m%d_%H%M%S"))

    if dpi is None:
        dpi = 300

    os.makedirs(path, exist_ok=True) # make dir if it doesn't already exist

    p9.save_as_pdf_pages(plots=[plot], filename=(filename), dpi=dpi)

    if verbose:
        print("Saved: {}".format(filename))