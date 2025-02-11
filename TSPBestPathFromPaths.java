import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPBestPathFromPaths {

    public static void main(String[] args) throws FileNotFoundException {
        File inputFile = new File("");
        Scanner scanner = new Scanner(inputFile);
        String[] input = scanner.nextLine().split(" ");
        double constant;
        int[][] distanceMatrix;
        int count = 0;
        if (input.length == 1) {
            constant = Double.parseDouble(input[0]);
            input = scanner.nextLine().split(" ");
            if (input[0].charAt(0) == '(') {
                ArrayList<String> inputData = new ArrayList<>();
                inputData.add(input[0] + " " + input[1]);
                while (scanner.hasNext()) {
                    inputData.add(scanner.nextLine().replace("\n", ""));
                }
                distanceMatrix = buildDistanceMatrixFromEdges(inputData, constant);
                count = distanceMatrix.length;
            } else {
                ArrayList<double[]> inputData = new ArrayList<>();
                double[] coordinates1 = new double[2];
                String[] temp;
                coordinates1[0] = Double.parseDouble(input[0]);
                coordinates1[1] = Double.parseDouble(input[1]);
                inputData.add(coordinates1);
                while (scanner.hasNext()) {
                    temp = scanner.nextLine().replace("\n", "").split(" ");
                    double[] coordinates = new double[2];
                    coordinates[0] = Double.parseDouble(temp[0]);
                    coordinates[1] = Double.parseDouble(temp[1]);
                    inputData.add(coordinates);
                }
                distanceMatrix = buildDistanceMatrixFromCoordinates(inputData, constant);
                count = distanceMatrix.length;
            }
        } else {
            List<int[]> edges = new ArrayList<>();
            edges.add(new int[]{Integer.parseInt(input[0]), Integer.parseInt(input[1]), (int) (Double.parseDouble(input[2]) * 100)});
            count = 1;
            while (scanner.hasNext()) {
                int city1 = scanner.nextInt();
                int city2 = scanner.nextInt();
                int cost = (int) (scanner.nextDouble() * 100);
                count++;
                edges.add(new int[]{city1, city2, cost});
            }
            int discriminant = 1 + 8 * count;
            int n1 = (int) ((1 + Math.sqrt(discriminant)) / 2);
            int n2 = (int) ((1 - Math.sqrt(discriminant)) / 2);
            count = Math.max(n1, n2);
            distanceMatrix = buildDistanceMatrixFromEdges(edges, count);
        }

        TSPBestPathFromPaths algorithm = new TSPBestPathFromPaths();
        List<Integer> bestPath = algorithm.FindTheBestPath(distanceMatrix);
        int bestCost = algorithm.calculatePathCost(bestPath, distanceMatrix);

        System.out.println("Best path found:");
        for (int city : bestPath) {
            System.out.print((city + 1) + " ");
        }
        System.out.println("\nBest cost: " + bestCost / 100.0);
    }

    public static int[][] buildDistanceMatrixFromEdges(List<int[]> edges, int numberOfCities) {
        int[][] distanceMatrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(distanceMatrix[i], Integer.MAX_VALUE / 2);
        }

        for (int[] edge : edges) {
            int city1 = edge[0];
            int city2 = edge[1];
            int cost = edge[2];
            distanceMatrix[city1][city2] = cost;
            distanceMatrix[city2][city1] = cost;
        }

        return distanceMatrix;
    }

    public List<Integer> FindTheBestPath(int[][] distanceMatrix) {
        int numCities = distanceMatrix.length;
        double[][] agentsStrength = new double[numCities][numCities];
        double[][] transitionProbabilities = new double[numCities][numCities];
        initializeAgentStrengths(agentsStrength);

        for (int i = 0; i < numCities; i++) {
            for (int j = 0; j < numCities; j++) {
                if (distanceMatrix[i][j] != 0) {
                    transitionProbabilities[i][j] = 1.0 / distanceMatrix[i][j];
                } else {
                    transitionProbabilities[i][j] = 0;
                }
            }
        }

        List<Integer> bestPath = null;
        int bestCost = Integer.MAX_VALUE;

        for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
            List<List<Integer>> allAgentsPaths = new ArrayList<>();
            List<Integer> allAgentsCosts = new ArrayList<>();

            for (int agent = 0; agent < NUM_AGENTS; agent++) {
                List<Integer> path = constructPath(distanceMatrix, agentsStrength, transitionProbabilities);
                int cost = calculatePathCost(path, distanceMatrix);
                allAgentsPaths.add(path);
                allAgentsCosts.add(cost);

                if (cost < bestCost) {
                    bestCost = cost;
                    bestPath = new ArrayList<>(path);
                }
            }

            updateAgentStrengths(agentsStrength, allAgentsPaths, allAgentsCosts);
            evaporateAgentStrengths(agentsStrength);

            System.out.println("Iteration " + (iteration + 1) + " completed. Best cost: " + bestCost);
        }

        return bestPath;
    }

    private void initializeAgentStrengths(double[][] agentsStrength) {
        for (int i = 0; i < agentsStrength.length; i++) {
            for (int j = 0; j < agentsStrength[i].length; j++) {
                agentsStrength[i][j] = 1.0;
            }
        }
    }

    private List<Integer> constructPath(int[][] distanceMatrix, double[][] agentsStrength, double[][] transitionProbabilities) {
        int numCities = distanceMatrix.length;
        List<Integer> path = new ArrayList<>();
        boolean[] visited = new boolean[numCities];

        int startCity = (int) (Math.random() * numCities);
        path.add(startCity);
        visited[startCity] = true;

        for (int i = 1; i < numCities; i++) {
            int currentCity = path.get(path.size() - 1);
            int nextCity = chooseNextCity(currentCity, visited, agentsStrength, transitionProbabilities);
            path.add(nextCity);
            visited[nextCity] = true;
        }

        return path;
    }

    private int chooseNextCity(int currentCity, boolean[] visited, double[][] agentsStrength, double[][] transitionProbabilities) {
        int numCities = agentsStrength.length;
        double[] probabilities = new double[numCities];
        double sum = 0.0;

        for (int city = 0; city < numCities; city++) {
            if (!visited[city]) {
                probabilities[city] = Math.pow(agentsStrength[currentCity][city], PARAM_A) *
                        Math.pow(transitionProbabilities[currentCity][city], PARAM_B);
                sum += probabilities[city];
            } else {
                probabilities[city] = 0;
            }
        }

        for (int city = 0; city < numCities; city++) {
            probabilities[city] /= sum;
        }

        double rand = Math.random();
        double cumulativeProbability = 0.0;
        for (int city = 0; city < numCities; city++) {
            cumulativeProbability += probabilities[city];
            if (rand <= cumulativeProbability) {
                return city;
            }
        }

        return -1;
    }

    private void updateAgentStrengths(double[][] agentsStrength, List<List<Integer>> allAgentsPaths, List<Integer> allAgentsCosts) {
        for (int i = 0; i < NUM_AGENTS; i++) {
            List<Integer> path = allAgentsPaths.get(i);
            int cost = allAgentsCosts.get(i);
            for (int j = 0; j < path.size() - 1; j++) {
                int cityA = path.get(j);
                int cityB = path.get(j + 1);
                agentsStrength[cityA][cityB] += FACTOR_D / cost;
                agentsStrength[cityB][cityA] += FACTOR_D / cost;
            }
        }
    }

    private void evaporateAgentStrengths(double[][] agentsStrength) {
        for (int i = 0; i < agentsStrength.length; i++) {
            for (int j = 0; j < agentsStrength[i].length; j++) {
                agentsStrength[i][j] *= (1 - PARAM_C);
            }
        }
    }

    public int calculatePathCost(List<Integer> path, int[][] distanceMatrix) {
        int cost = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            cost += distanceMatrix[path.get(i)][path.get(i + 1)];
        }
        cost += distanceMatrix[path.get(path.size() - 1)][path.get(0)];
        return cost;
    }

    public static int calculateDistance(double x1, double y1, double x2, double y2) {
        return (int) ((Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2))) * 10);
    }

    public static int[][] buildDistanceMatrixFromCoordinates(List<double[]> coordinates, double cost) {
        int numberOfCities = coordinates.size();
        int[][] distanceMatrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(distanceMatrix[i], Integer.MAX_VALUE);
        }

        for (int i = 0; i < numberOfCities; i++) {
            for (int j = 0; j < numberOfCities; j++) {
                if (i != j) {
                    double[] point1 = coordinates.get(i);
                    double[] point2 = coordinates.get(j);
                    double distance = calculateDistance(point1[0], point1[1], point2[0], point2[1]);
                    distanceMatrix[i][j] = (int) (distance * cost);
                    distanceMatrix[j][i] = (int) (distance * cost);
                } else {
                    distanceMatrix[i][j] = 0;
                }
            }
        }
        return distanceMatrix;
    }

    public static int[][] buildDistanceMatrixFromEdges(List<String> edgesInput, double constant) {
        Map<Integer, double[]> cityCoordinates = new HashMap<>();
        Map<String, Integer> coordinateToIndex = new HashMap<>();

        int cityIndex = 0;

        List<double[]> edges = new ArrayList<>();
        for (String edgeInput : edgesInput) {
            String[] points = edgeInput.split(" ");
            String point1Str = points[0].replace("(", "").replace(")", "");
            String point2Str = points[1].replace("(", "").replace(")", "");

            double[] point1 = Arrays.stream(point1Str.split(",")).mapToDouble(Double::parseDouble).toArray();
            double[] point2 = Arrays.stream(point2Str.split(",")).mapToDouble(Double::parseDouble).toArray();

            if (!coordinateToIndex.containsKey(point1Str)) {
                coordinateToIndex.put(point1Str, cityIndex);
                cityCoordinates.put(cityIndex, point1);
                cityIndex++;
            }
            if (!coordinateToIndex.containsKey(point2Str)) {
                coordinateToIndex.put(point2Str, cityIndex);
                cityCoordinates.put(cityIndex, point2);
                cityIndex++;
            }

            int city1 = coordinateToIndex.get(point1Str);
            int city2 = coordinateToIndex.get(point2Str);

            double distance = calculateDistance(point1[0], point1[1], point2[0], point2[1]);

            edges.add(new double[]{city1, city2, distance});
        }

        int numberOfCities = cityIndex;
        int[][] distanceMatrix = new int[numberOfCities][numberOfCities];
        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(distanceMatrix[i], Integer.MAX_VALUE);
        }

        for (double[] edge : edges) {
            int city1 = (int) edge[0];
            int city2 = (int) edge[1];
            double distance = edge[2];
            distanceMatrix[city1][city2] = (int) (constant * distance);
            distanceMatrix[city2][city1] = (int) (constant * distance);
        }

        return distanceMatrix;
    }

    private static final int PARAM_A = 1;
    private static final int PARAM_B = 5;
    private static final double PARAM_C = 0.5;
    private static final int FACTOR_D = 100;
    private static final int NUM_AGENTS = 20;
    private static final int MAX_ITERATIONS = 1000;
}
