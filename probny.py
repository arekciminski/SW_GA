# importing libraries
import numpy as np
import time
import matplotlib.pyplot as plt

# creating initial data values
# of x and y
x = np.linspace(0, 10, 100)
y = np.sin(x)
y1 = np.sin(2*x)

# to run GUI event loop
plt.ion()

# here we are creating sub plots
figure, ax = plt.subplots(figsize=(10, 8))
line1, = ax.plot(x, y)
line2, = ax.plot(x, y1)

# setting title
plt.title("Geeks For Geeks", fontsize=20)

# setting x-axis label and y-axis label
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

# Loop
for i in range(50):
    # creating new Y values
    new_y = np.sin(x - 0.5 * i)
    new_y2 = i*np.sin(2*x - 0.5 * i)

    # updating data values
    line1.set_xdata(x)
    line1.set_ydata(new_y)
    line2.set_xdata(x)
    line2.set_ydata(new_y2)

    # drawing updated values
    figure.canvas.draw()

    # This will run the GUI event
    # loop until all UI events
    # currently waiting have been processed
    figure.canvas.flush_events()

    time.sleep(0.1)