#!/usr/local/bin/python
'''

 Simple illustration for plotting coast lines in python


'''

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import scipy.io          as spio
import matplotlib.pyplot as plt

################################################################################
#                    E X P O R T E D     C L A S S E S:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#             U N I T     T E S T     C A S E     F U N C T I O N:
################################################################################

#-------------------------------------------------------------------------------


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
  '''
  Simple illustration for plotting coast lines in python

  '''

  #  Set the name of the file with the coastline data
  coastFile = 'earth_coastline.mat'

  #  Load the file.
  earth_coastline = spio.loadmat( coastFile )['earth_coastline']

  #  Generate the figure for the plot
  plt.figure()

  #  Plot the coastlines
  plt.plot( earth_coastline[:,0], earth_coastline[:,1], 'k' )

  #  Set the correct aspect ratio for the plot
  plt.axis( 'square' )
  plt.axis( [ -180, 180, -90, 90 ] )

  #  Show the plot
  plt.show()

  return

if __name__ == "__main__":
  main()
