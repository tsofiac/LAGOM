import pandas as pd
import matplotlib.pyplot as plt

def plot_and_save_loss(csv_file_path, output_file_path):
    # Read the CSV file
    data = pd.read_csv(csv_file_path, sep='\t')  # Use tab separator for your data

    epoch = data['epoch']
    training_loss = data['training_loss']
    validation_loss = data['validation_loss']

    # Plot the training and validation loss against epochs
    plt.figure(figsize=(10, 6))
    
    plt.plot(epoch, training_loss, label='Training Loss', color='blue')
    plt.plot(epoch, validation_loss, label='Validation Loss', color='orange')
    
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss vs. Epochs')
    plt.legend()
    plt.grid(True)
    
    # Show the plot
    plt.show()

 
    # Save the figure
    plt.savefig(output_file_path, format='png')

    # Show the plot
    plt.show()

# Example usage
csv_file_path = 'Sofia_and_Miranda/results/training/backward_prediction/version_1/logged_train_metrics.csv'
output_file_path = 'loss_plot.png'
plot_and_save_loss(csv_file_path, output_file_path)