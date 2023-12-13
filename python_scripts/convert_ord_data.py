import tensorflow as tf
import csv

def convert_ord_data_to_csv(model_path):
    # Load the SavedModel
    model = tf.saved_model.load(model_path)

    # Create a CSV file for writing
    with open('output.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write header row with variable names
        header = ["VariableName", "Value"]
        csv_writer.writerow(header)

        # Iterate through the model's variables
        for variable in model.variables:
            # Get the variable name
            variable_name = variable.name

            # Get the variable value as a NumPy array
            variable_value = variable.numpy()

            # Write the variable name and value to the CSV file
            csv_writer.writerow([variable_name, variable_value])

# Specify the path to the SavedModel
model_path = '/Users/ashort/project_1/python_scripts/0a/SavedModel'

# Call the function with the path
convert_ord_data_to_csv(model_path)




