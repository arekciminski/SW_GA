'''# importing libraries
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

    time.sleep(0.1)'''
# SuperFastPython.com
# example of a parallel for loop with multiple arguments
from time import sleep
from random import random
from multiprocessing import Pool


# task to execute in another process
def task(arg1, arg2, arg3):
    # generate a value between 0 and 1
    value = random()
    # block for a fraction of a second to simulate work
    sleep(value)
    # return the generated value
    return value


# entry point for the program
if __name__ == '__main__':
    # create the process pool
    with Pool() as pool:
        # prepare arguments
        items = [[i, i * 2, i * 3] for i in range(10)]
        print(items)
        # call the same function with different data in parallel
        for result in pool.starmap(task, items):
            # report the value to show progress
            print(result)