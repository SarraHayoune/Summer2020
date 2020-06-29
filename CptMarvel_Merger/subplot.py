import matplotlib.pyplot as plt
import matplotlib.image as mpimg


img1= mpimg.imread('ID= 89425759,z= 3.0.png')
img2= mpimg.imread('ID= 89425759,z= 2.65.png')
img3= mpimg.imread('ID= 89425759,z= 2.29.png')

img4= mpimg.imread('ID= 89425759,z= 1.75.png')
img5= mpimg.imread('ID= 89425759,z= 1.5.png')
img6= mpimg.imread('ID= 89425759,z= 1.4.png')

#(ax4,ax5,ax6)

fig,((ax1, ax2,ax3))= plt.subplots(1,3, figsize=(10,10))
ax1.imshow(img1)
ax1.axis('off')
ax2.imshow(img2)
ax2.axis('off')
ax3.imshow(img3)
ax3.axis('off')
#ax4.imshow(img4)
#ax4.axis('off')
#ax5.imshow(img5)
#ax5.axis('off')
#ax6.imshow(img6)
#ax6.axis('off')
fig.tight_layout()

fig.savefig('merger1.png')
