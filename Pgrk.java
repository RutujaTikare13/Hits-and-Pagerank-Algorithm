// RUTUJA_6153
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class Pgrk6153 {
    private static final double dampingFactor = 0.85;
    double staticScore = 0.0;
    int iterations;
    double initialValue;
    int vertexCount;
    int edgeCount;
    double errorRate;
    int count = 0;
    boolean isErrorRateTermination = true;
    boolean isTerminated = false;
    boolean isLargeGraph = false;
    double[][] adjacencyMatrix;
    double[][] adjacencyMatrixTranspose;
    double[] pageRantVector;


    public static void main(String[] args) {
        Pgrk6153 pgrk6153 = new Pgrk6153();

        if (args.length != 3) {
            System.out.println("Invalid number of arguments. The correct format is:" +
                    "java Pgrk6153 <iterations> <initialValue> <fileName>");
            return;
        }
        pgrk6153.iterations = Integer.parseInt(args[0]);
        pgrk6153.initialValue = Double.parseDouble(args[1]);

        try {
            pgrk6153.init(args[2]);
        } catch (Exception e) {
            System.out.println("Error in execution: " + e.getMessage());
            return;
        }
        pgrk6153.initializeVectors();
        if (!pgrk6153.isLargeGraph) {
            pgrk6153.displayUpdatedVectors();
        }

        while (!pgrk6153.isTerminated) {
            pgrk6153.updatePageRank();
        }

        if (pgrk6153.isLargeGraph) {
            pgrk6153.displayLastIteration();
        }
    }


    private void checkTermination(double[] tempPageRank) {
        if (!isErrorRateTermination) {
            if(count >= iterations) {
                isTerminated = true;
            }
        } else {
            isTerminated = true;
            for(int i = 0; i < tempPageRank.length; i++) {
                if (Math.abs(tempPageRank[i] - pageRantVector[i]) >= errorRate) {
                    isTerminated = false;
                    break;
                }
            }
        }
    }

    private void updatePageRank() {
        double[] tempPageRankVector = matrixVectorDotProduct(adjacencyMatrixTranspose, pageRantVector);

        count++;
        checkTermination(tempPageRankVector);
        pageRantVector = tempPageRankVector;
        if (!isLargeGraph) {
            displayUpdatedVectors();
        }
    }

    private void displayLastIteration() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("Iter")
                .append("\t")
                .append(": ")
                .append(count)
                .append("\n");
        for (int i = 0; i < pageRantVector.length; i++) {
            stringBuilder.append("A/H[")
                    .append(i).append("]=")
                    .append(String.format("%.7f", pageRantVector[i]))
                    .append("\n");
        }
        System.out.println(stringBuilder.toString());
    }

    private void displayUpdatedVectors() {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(count == 0 ? "Base" : "Iter")
                .append("\t")
                .append(count)
                .append(" :");
        for (int i = 0; i < pageRantVector.length; i++) {
            stringBuilder.append("P[")
                    .append(i)
                    .append("]=")
                    .append(String.format("%.7f", pageRantVector[i]))
                    .append(" ");
        }
        System.out.println(stringBuilder.toString());

    }

    private double[] matrixVectorDotProduct(double[][] matrix, double[] vector) {
        double[] updatedVector = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            double currentScore = 0.0;
            for (int j = 0; j < vector.length; j++) {
                currentScore += (matrix[i][j] * vector[j]);
            }
            updatedVector[i] = staticScore + (dampingFactor * currentScore);
        }
        return updatedVector;
    }

    private void initializeVectors() {
        Arrays.fill(pageRantVector, initialValue);
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
                    adjacencyMatrix = new double[left][left];
                    adjacencyMatrixTranspose = new double[left][left];
                    pageRantVector = new double[left];
                    for(int i = 0; i < adjacencyMatrix.length; i++) {
                        Arrays.fill(adjacencyMatrix[i], 0);
                    }
                    vertexCount = left;
                    edgeCount = right;
                    isFirstLine = false;
                    continue;
                }
                adjacencyMatrix[left][right] = 1.0;
                recordCount++;
            }

            for (int i = 0; i < vertexCount; i++) {
                int counter = 0;
                for (int j = 0; j < vertexCount; j++) {
                    if (adjacencyMatrix[i][j] == 1) {
                        counter++;
                    }
                }
                if (counter > 0) {
                    for (int j = 0; j < vertexCount; j++) {
                            adjacencyMatrix[i][j] = (adjacencyMatrix[i][j])/counter;
                    }
                }
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
                throw new RuntimeException("Number of Edged Mismatch");
            }

            staticScore = (1.0 - dampingFactor)/ vertexCount;

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
