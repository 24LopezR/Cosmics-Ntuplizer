import ROOT as R
def setStyle(width=800, height=600, font=42, fontsize=0.04):
    style = R.gStyle

    width = width
    height = height
    font = font
    tMargin = 0.1 
    lMargin = 0.15
    fontsize = float(fontsize)

    rMargin = tMargin * float(height) / float(width)
    bMargin = lMargin
    titleX = lMargin + (1-lMargin-rMargin)/2
    titleY = 1 - (tMargin/2)

    # canvas
    style.SetCanvasDefW(width)               # width
    style.SetCanvasDefH(height)              # height

    # pad margins
    style.SetPadTopMargin(tMargin)           # default 0.1
    style.SetPadBottomMargin(bMargin)        # default 0.1
    style.SetPadLeftMargin(lMargin)          # default 0.1
    style.SetPadRightMargin(rMargin)         # default 0.1

    # legend
    style.SetLegendFont(font)                # helvetica normal
    style.SetLegendTextSize(fontsize)        # default 0
    style.SetLegendBorderSize(0)             # off

    # title
    style.SetTitleFont(font, '')              # helvetica normal
    style.SetTitleFontSize(fontsize)         # default 0
    style.SetTitleX(titleX)                  # center title horizontally with respect to frame
    style.SetTitleY(titleY)                  # center title vertically within margin

    # axis titles
    style.SetTitleFont(font, 'XYZ')          # helvetica normal
    style.SetTitleSize(fontsize, 'XYZ')      # default 0.02

    # axis labels
    style.SetLabelFont(font, 'XYZ')          # helvetica normal
    style.SetLabelSize(fontsize, 'XYZ')      # default 0.04

    style.cd()
