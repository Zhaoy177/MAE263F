import numpy as np
import matplotlib.pyplot as plt

# Seed for reproducibility
np.random.seed(42)

# Generate 10 random x values within a range
x_generated = np.linspace(0, 5, 10)

# Parameters for the function (can use the previously fitted values or set randomly)
n_true = 0.06
a_true = 0.25
m_true = 0.57
b_true = 0.11

# Generate corresponding y values based on the function with added noise
noise = 0.001 * np.random.normal(0, 0.1, size=x_generated.shape)  # Add Gaussian noise
y_generated = n_true * np.exp(-a_true * (m_true * x_generated + b_true) ** 2) + noise

# Display the generated x and y arrays
x_data = x_generated
y_data = y_generated

def compute_loss(x, y, m, b):
    """
    Compute the Mean Squared Error (MSE) loss for the linear function y = mx + b.

    Parameters:
    x : np.array
        Input data points (x values).
    y : np.array
        Actual output data points (y values).
    m, b : float
        Parameters of the linear function.

    Returns:
    float
        Mean Squared Error (MSE) loss.
    """
    y_pred = m * x + b
    return np.mean((y - y_pred) ** 2)

# Hyperparameter search
learning_rates = np.arange(0.0001, 0.01, 0.0002)
epochs_list = np.arange(1000, 100001, 2000)

best_loss = float('inf')
best_params = None

for learning_rate in learning_rates:
    for epochs in epochs_list:
        # Reinitialize parameters
        m_fit = np.random.rand()
        b_fit = np.random.rand()

        for epoch in range(epochs):
            # Forward pass: compute predicted outputs
            y_pred_fit = m_fit * x_data + b_fit

            # Compute gradients
            grad_m_fit = -2 * np.mean((y_data - y_pred_fit) * x_data)
            grad_b_fit = -2 * np.mean(y_data - y_pred_fit)

            # Update parameters
            m_fit -= learning_rate * grad_m_fit
            b_fit -= learning_rate * grad_b_fit

        # Compute final loss
        loss_fit = compute_loss(x_data, y_data, m_fit, b_fit)

        # Save best parameters
        if loss_fit < best_loss:
            best_loss = loss_fit
            best_params = (learning_rate, epochs, m_fit, b_fit)

        print(f"Learning Rate: {learning_rate}, Epochs: {epochs}, Loss: {loss_fit:.6f}")

# Display the best parameters and their corresponding loss
print("\nBest Combination:")
print(f"Learning Rate: {best_params[0]}, Epochs: {best_params[1]}")
print(f"m: {best_params[2]:.6f}, b: {best_params[3]:.6f}")
print(f"Loss: {best_loss:.6f}")

# Predicted y values using the best parameters
y_predicted = best_params[2] * x_generated + best_params[3]

# Plot the training data
plt.scatter(x_generated, y_generated, color='blue', label='Training Data (Noisy)', marker='o')

# Plot the predicted data
plt.plot(x_generated, y_predicted, color='red', label='Predicted Data (Model)', linestyle='--')

# Add labels, title, and legend
plt.xlabel("x")
plt.ylabel("y")
plt.title("Comparison of Training Data and Predicted Data")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
