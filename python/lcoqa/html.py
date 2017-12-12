import os
import astropy.time as time
import sdss_access

spath = sdss_access.path.Path()


class Summary(object):
    """Object for writing HTML summary web pge

    Attributes:
    ----------
    mjd : np.int32, int
      integer MJD
    htmlfile : 
    body_text : 
    header_text : 
    framework : 

    Methods:
    -------
    add_framework : 
    add_exposure : 
    add_exposure_table : 
    write : 
    _cell_plot : 

    """
    def __init__(self, mjd=None, htmlfile=None):
        self.header_text = ""
        if(mjd is not None):
            tt = time.Time(mjd, format='mjd', out_subfmt='date')
            tt.format = 'isot'
            self.header_text = """<h1>{mjd}<br/>
{date}
</h1>
""".format(date=tt, mjd=mjd)
        self.body_text = ""
        self.mjd = mjd
        self.htmlfile = htmlfile
        self.add_framework()
        self.add_exposure_table()
        return

    def add_framework(self):
        self.framework = """<html>
  <header>
     {header_text}
  </header>
  <body>
     {body_text}
  </body>
</html>
"""
        return

    def add_exposure_table(self):
        self.body_text += """
    <table>
{exposure_table_text}
    </table>
"""
        self.exposure_table_text = ""
        return

    def _cell_plot(self, prefix=None, plotdir=None, mjd=None, plate=None,
                   expno=None):
        pngfile = "{plotdir}/{prefix}-{mjd}-{plate}-{expno}.png'"
        pngfile = pngfile.format(prefix=prefix, plotdir=plotdir,
                                 mjd=mjd, plate=plate, expno=expno)
        cell = "<a href='{pngfile}'><img src='{pngfile}' width=300px /></a>"
        cell = cell.format(pngfile=pngfile)
        return cell

    def add_exposure(self, exposure=None):
        plotdir = os.path.join(str(exposure['mjd']), 'plots')
        row = "<tr>"
        os.environ['PLATELIST_DIR'] = "http://platedesign.sdss.org"
        plink = os.path.dirname(spath.full('plateHoles',
                                           plateid=exposure['plateid']))
        row += "<td><a href={plink}><b>PLATE={plate}</b></a><br/>".format(plate=exposure['plateid'], plink=plink)
        row += "EXPOSURE={expno}<br/>".format(expno=exposure['expno'])
        row += "RA={raCen:>0.3f}<br/>".format(raCen=exposure['raCen'])
        row += "DEC={decCen:>0.3f}<br/>".format(decCen=exposure['decCen'])
        row += "MJD={mjd}<br/>".format(mjd=exposure['mjd'])
        row += "CART={cart}<br/>".format(cart=exposure['cartid'])
        row += "DATE={dateobs}<br/>".format(dateobs=exposure['date-obs'])
        row += "HA={ha:>0.2f}<br/>".format(ha=exposure['ha'])
        row += "HA_DESIGN={ha_design:>0.2f}<br/>".format(ha_design=exposure['ha_design'])
        row += "HA_MIN={ha_observable_min:>0.2f}<br/>".format(ha_observable_min=exposure['ha_observable_min'])
        row += "HA_MAX={ha_observable_max:>0.2f}<br/>".format(ha_observable_max=exposure['ha_observable_max'])
        row += "SEEING={seeing:>0.2f}<br/>".format(seeing=exposure['seeing_median'])
        row += "GUIDERMS={guiderms:>0.2f}</td>".format(guiderms=exposure['gdrms_median'])
        seeing_cell = self._cell_plot(prefix='seeing', plotdir=plotdir,
                                      mjd=exposure['mjd'],
                                      plate=exposure['plateid'],
                                      expno=exposure['expno'])
        gdrms_cell = self._cell_plot(prefix='guider-rms', plotdir=plotdir,
                                     mjd=exposure['mjd'],
                                     plate=exposure['plateid'],
                                     expno=exposure['expno'])
        signal_cell = self._cell_plot(prefix='signal', plotdir=plotdir,
                                      mjd=exposure['mjd'],
                                      plate=exposure['plateid'],
                                      expno=exposure['expno'])
        sn2_cell = self._cell_plot(prefix='sn2', plotdir=plotdir,
                                   mjd=exposure['mjd'],
                                   plate=exposure['plateid'],
                                   expno=exposure['expno'])
        signal_focal_cell = self._cell_plot(prefix='signal-focal',
                                            plotdir=plotdir,
                                            mjd=exposure['mjd'],
                                            plate=exposure['plateid'],
                                            expno=exposure['expno'])
        sn2_focal_cell = self._cell_plot(prefix='sn2-focal',
                                         plotdir=plotdir,
                                         mjd=exposure['mjd'],
                                         plate=exposure['plateid'],
                                         expno=exposure['expno'])
        row += "<td>{signal_cell}</td>".format(signal_cell=signal_cell)
        row += "<td>{sn2_cell}</td>".format(sn2_cell=sn2_cell)
        row += "<td>{signal_focal_cell}</td>".format(signal_focal_cell=signal_focal_cell)
        row += "<td>{sn2_focal_cell}</td>".format(sn2_focal_cell=sn2_focal_cell)
        row += "<td>{seeing_cell}</td>".format(seeing_cell=seeing_cell)
        row += "<td>{gdrms_cell}</td>".format(gdrms_cell=gdrms_cell)
        row += "</tr>\n"
        self.exposure_table_text += row

    def write(self):
        body_text_final = self.body_text.format(exposure_table_text=self.exposure_table_text)
        text = self.framework.format(header_text=self.header_text,
                                     body_text=body_text_final)
        fp = open(self.htmlfile, "w")
        fp.write(text)
        fp.close()
