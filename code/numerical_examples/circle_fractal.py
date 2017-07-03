from math import sin, cos, pi

shift = pi/5
def getcirc(x, y, r, level, base=pi/2):
    cs = [[x, y, r]]
    if level > 0:
        cs += getcirc(x + r*cos(base+shift), y + r*sin(base+shift), r/1.5, level-1, base+shift)
        cs += getcirc(x + r*cos(base-shift), y + r*sin(base-shift), r/1.5, level-1, base-shift)
    return cs

crs = getcirc(0, 0, 0.5, 4)

#  import matplotlib.pyplot as plt

#  for c in crs:
#      circle=plt.Circle((c[0],c[1]),c[2],color='r')
#      plt.gcf().gca().add_artist(circle)
#
#  plt.show()

print(str(crs).replace('[','{').replace(']','}'))
