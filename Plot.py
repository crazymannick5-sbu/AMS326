import matplotlib.pyplot as plt
import csv

def plot_csv(file_name):
    x_values = []
    y_values = []

    with open(file_name, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            x_values.append(float(row[0]))
            y_values.append(float(row[1]))

    plt.plot(x_values, y_values)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Data Plot')
    plt.show()

# Example usage:
file_name = 'Data.csv'
plot_csv(file_name)
