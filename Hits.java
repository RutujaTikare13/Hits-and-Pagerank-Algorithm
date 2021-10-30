// RUTUJA_6153
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class Hits6153 {
    int iterations;
    double initialValue;
    int vertexCount;
    int edgeCount;
    double errorRate;
    int count = 0;
    boolean isErrorRateTermination = true;
    boolean isTerminated = false;
    boolean isLargeGraph = false;
    int[][] adjacencyMatrix;
    int[][] adjacencyMatrixTranspose;
    double[] hubVector;
    double[] authVector;

    public static void main(String[] args) {
        Hits6153 hits6153 = new Hits6153();
        if (args.length != 3) {
            System.out.println("Invalid number of arguments. The correct format is:" +
                    "java Hits6153 <iterations> <initialValue> <fileName>");
            return;
        }
        try {
            hits6153.iterations = Integer.parseInt(args[0]);
            hits6153.initialValue = Double.parseDouble(args[1]);
        } catch (Exception e) {
            System.out.println("Invalid iterations or initial value");
            return;
        }

        try {
            hits6153.init(args[2]);
        } catch (Exception e) {
            System.out.println("Error in execution: " + e.getMessage());
            return;
        }
        hits6153.initializeVectors();
        if (!hits6153.isLargeGraph) {
            hits6153.displayUpdatedVectors();
        }
        while (!hits6153.isTerminated) {
            hits6153.updateHubAndAuth();
        }
        if (hits6153.isLargeGraph) {
            hits6153.displayLastIteration();
        }
    }

    private void checkTermination(double[] tempAuthVector, double[] tempHubVector) {
        if (!isErrorRateTermination) {
            if(count >= iterations) {
                isTerminated = true;
            }
        } else {
            isTerminated = true;
            for(int i = 0; i < tempAuthVector.length; i++) {
                if (Math.abs(tempAuthVector[i] - authVector[i]) >= errorRate ||
                        Math.abs(tempHubVector[i] - hubVector[i]) >= errorRate) {
                    isTerminated = false;
                    break;
                }
            }
        }
    }

    private void updateHubAndAuth() {
        double[] tempAuthVector = matrixVectorDotProduct(adjacencyMatrixTranspose, hubVector);
        double[] tempHubVector = matrixVectorDotProduct(adjacencyMatrix, tempAuthVector);
        double rootMeanSquareHub = 0;
        double rootMeanSquareAuth = 0;

        for(int i = 0; i < tempAuthVector.length; i ++) {
            rootMeanSquareHub += (tempHubVector[i] * tempHubVector[i]);
            rootMeanSquareAuth += (tempAuthVector[i] * tempAuthVector[i]);
        }
        rootMeanSquareHub = Math.sqrt(rootMeanSquareHub);
        rootMeanSquareAuth = Math.sqrt(rootMeanSquareAuth);

        for(int i = 0; i < tempAuthVector.length; i++) {
            tempHubVector[i] = tempHubVector[i] / rootMeanSquareHub;
            tempAuthVector[i] = tempAuthVector[i] / rootMeanSquareAuth;
        }
        count++;
        checkTermination(tempAuthVector, tempHubVector);
        authVector = tempAuthVector;
        hubVector = tempHubVector;
        if (!isLargeGraph) {
            displayUpdatedVectors();
        }
    }

    private void displayLastIteration() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("Iter").append("\t").append(": ").append(count).append("\n");
        for (int i = 0; i < authVector.length; i++) {
            stringBuilder.append("A/H[")
                    .append(i)
                    .append("]=")
                    .append(String.format("%.7f", authVector[i]))
                    .append("/")
                    .append(String.format("%.7f", hubVector[i]))
                    .append("\n");
        }
        System.out.println(stringBuilder.toString());
    }

    private void displayUpdatedVectors() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(count == 0 ? "Base" : "Iter")
                .append("\t")
                .append(count).append(" :");
        for (int i = 0; i < authVector.length; i++) {
            stringBuilder.append("A/H[")
                    .append(i)
                    .append("]=")
                    .append(String.format("%.7f", authVector[i]))
                    .append("/")
                    .append(String.format("%.7f", hubVector[i]))
                    .append(" ");
        }
        System.out.println(stringBuilder.toString());
    }

    private double[] matrixVectorDotProduct(int[][] matrix, double[] vector) {
        double[] updatedVector = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < vector.length; j++) {
                updatedVector[i] += matrix[i][j] * vector[j];
            }
        }
        return updatedVector;
    }

    private void initializeVectors() {
        Arrays.fill(hubVector, initialValue);
        Arrays.fill(authVector, initialValue);
    }

    private void init(String fileName) {
        BufferedReader reader;
        String input;
        int left;
        int right;
        boolean isFirstLine = true;
        int recordCount = 0;
        try {
            reader = new BufferedReader(new FileReader(fileName));
            while ((input = reader.readLine()) != null) {
                String[] splits = input.split(" ");
                left = Integer.parseInt(splits[0]);
                right = Integer.parseInt(splits[1]);
                if (isFirstLine) {
                    adjacencyMatrix = new int[left][left];
                    adjacencyMatrixTranspose = new int[left][left];
                    hubVector = new double[left];
                    authVector = new double[left];
                    for(int i = 0; i < adjacencyMatrix.length; i++) {
                        Arrays.fill(adjacencyMatrix[i], 0);
                    }
                    vertexCount = left;
                    edgeCount = right;
                    isFirstLine = false;
                    continue;
                }
                adjacencyMatrix[left][right] = 1;
                recordCount++;
            }
            for(int i = 0; i < vertexCount; i++) {
                for (int j = 0; j < vertexCount; j++) {
                    adjacencyMatrixTranspose[i][j] = adjacencyMatrix[j][i];
                }
            }
            if (vertexCount > 10) {
                iterations = 0;
                initialValue = -1;
                isLargeGraph = true;
            }

            if(recordCount != edgeCount) {
                throw new RuntimeException("Number of Edges Mismatch");
            }

            if (iterations == 0) {
                errorRate = Math.pow(10, -5);
                isErrorRateTermination = true;
            } else if (iterations < 0 && iterations >= -6) {
                errorRate = Math.pow(10, iterations);
                isErrorRateTermination = true;
            } else if (iterations < -6) {
                throw new RuntimeException("Iteration can not be below -6");
            } else {
                isErrorRateTermination = false;
            }

            if (!(initialValue == 1.0 || initialValue == 0.0 || initialValue == -1.0 || initialValue == -2.0)) {
                throw new RuntimeException("Invalid Initial Value. Initial value can be -2 or -1 or 0 or 1");
            }

            if (initialValue == -1) {
                initialValue = 1.0/ vertexCount;
            } else if (initialValue == -2) {
                initialValue = 1.0 / Math.sqrt(vertexCount);
            }
        } catch (FileNotFoundException e) {
            System.out.println("File: " + fileName + " Not Found");
            throw new RuntimeException(e);
        } catch (IOException e) {
            System.out.println("Unable to read contents of the file: " + fileName + ". Exiting!!!.");
            throw new RuntimeException(e);
        } catch (Exception e) {
            System.out.println("Exception while executing.");
            throw new RuntimeException(e);
        }
    }
}
