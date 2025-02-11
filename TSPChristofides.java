import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPChristofides {

    public static int[][] buildCostMatrix(List<int[]> edges, int numberOfCities) {
        int[][] costMatrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Integer.MAX_VALUE / 2);
        }

        for (int[] edge : edges) {
            int city1 = edge[0];
            int city2 = edge[1];
            int cost = edge[2];
            costMatrix[city1][city2] = cost;
            costMatrix[city2][city1] = cost;
        }

        return costMatrix;
    }
    public static int[][] buildCostMatrixFromCoordinates(List<double[]> coordinates, double cost) {
        int numberOfCities = coordinates.size();
        int[][] costMatrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Integer.MAX_VALUE);
        }

        for (int i = 0; i < numberOfCities; i++) {
            for (int j = 0; j < numberOfCities; j++) {
                if (i != j) {
                    double[] point1 = coordinates.get(i);
                    double[] point2 = coordinates.get(j);
                    double distance = calculateDistance(point1[0], point1[1], point2[0], point2[1]);
                    costMatrix[i][j] = (int) (distance * cost);
                    costMatrix[j][i] = (int) (distance * cost);
                } else {
                    costMatrix[i][j] = 0;
                }
            }
        }
        return costMatrix;
    }
    public static int[][] buildCostMatrix2(List<String> edgesInput, double constant) {
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
        int[][] costMatrix = new int[numberOfCities][numberOfCities];
        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Integer.MAX_VALUE);
        }

        for (double[] edge : edges) {
            int city1 = (int) edge[0];
            int city2 = (int) edge[1];
            double distance = edge[2];
            costMatrix[city1][city2] = (int) (constant * distance);
            costMatrix[city2][city1] = (int) (constant * distance);
        }

        return costMatrix;
    }

    private List<Integer> findOddDegreeVertices(List<int[]> mst, int numberOfCities) {
        int[] degree = new int[numberOfCities];
        for (int[] edge : mst) {
            degree[edge[0]]++;
            degree[edge[1]]++;
        }

        List<Integer> oddVertices = new ArrayList<>();
        for (int i = 0; i < numberOfCities; i++) {
            if (degree[i] % 2 != 0) {
                oddVertices.add(i);
            }
        }

        return oddVertices;
    }
    private List<int[]> findMinimumWeightMatching(List<Integer> oddVertices, int[][] costMatrix) {
        List<int[]> matching = new ArrayList<>();
        boolean[] matched = new boolean[costMatrix.length];

        for (int i = 0; i < oddVertices.size(); i++) {
            if (matched[oddVertices.get(i)]) continue;

            int bestMatch = -1;
            int minCost = Integer.MAX_VALUE;
            for (int j = i + 1; j < oddVertices.size(); j++) {
                if (!matched[oddVertices.get(j)] && costMatrix[oddVertices.get(i)][oddVertices.get(j)] < minCost) {
                    bestMatch = oddVertices.get(j);
                    minCost = costMatrix[oddVertices.get(i)][oddVertices.get(j)];
                }
            }

            if (bestMatch != -1) {
                matched[oddVertices.get(i)] = true;
                matched[bestMatch] = true;
                matching.add(new int[]{oddVertices.get(i), bestMatch, minCost});
            }
        }

        return matching;
    }
    private ArrayList<Integer> findEulerianCircuit(List<int[]> multigraph, int numberOfCities) {
        Map<Integer, List<Integer>> adjacencyList = new HashMap<>();
        for (int i = 0; i < numberOfCities; i++) {
            adjacencyList.put(i, new ArrayList<>());
        }
        for (int[] edge : multigraph) {
            adjacencyList.get(edge[0]).add(edge[1]);
            adjacencyList.get(edge[1]).add(edge[0]);
        }

        Stack<Integer> stack = new Stack<>();
        ArrayList<Integer> circuit = new ArrayList<>();
        stack.push(0);

        while (!stack.isEmpty()) {
            int current = stack.peek();
            if (adjacencyList.get(current).isEmpty()) {
                circuit.add(current);
                stack.pop();
            } else {
                int next = adjacencyList.get(current).remove(0);
                adjacencyList.get(next).remove((Integer) current);
                stack.push(next);
            }
        }

        return circuit;
    }
    private List<int[]> findMST(int[][] costMatrix) {
        int n = costMatrix.length;
        boolean[] inMST = new boolean[n];
        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(edge -> edge[2]));
        List<int[]> mst = new ArrayList<>();

        inMST[0] = true;
        for (int i = 1; i < n; i++) {
            if (costMatrix[0][i] < Integer.MAX_VALUE / 2) {
                pq.offer(new int[]{0, i, costMatrix[0][i]});
            }
        }

        while (!pq.isEmpty() && mst.size() < n - 1) {
            int[] edge = pq.poll();
            int u = edge[0];
            int v = edge[1];
            if (inMST[v]) continue;

            inMST[v] = true;
            mst.add(edge);

            for (int w = 0; w < n; w++) {
                if (!inMST[w] && costMatrix[v][w] < Integer.MAX_VALUE / 2) {
                    pq.offer(new int[]{v, w, costMatrix[v][w]});
                }
            }
        }

        return mst;
    }

    private ArrayList<Integer> convertToHamiltonianCircuit(List<Integer> eulerianCircuit) {
        Set<Integer> visited = new HashSet<>();
        ArrayList<Integer> hamiltonianCircuit = new ArrayList<>();

        for (int city : eulerianCircuit) {
            if (!visited.contains(city)) {
                visited.add(city);
                hamiltonianCircuit.add(city);
            }
        }

        hamiltonianCircuit.add(hamiltonianCircuit.get(0));
        return hamiltonianCircuit;
    }

    public int calculateTourCost(List<Integer> tour, int[][] costMatrix) {
        int cost = 0;
        for (int i = 0; i < tour.size() - 1; i++) {
            cost += costMatrix[tour.get(i)][tour.get(i + 1)];
        }
        cost += costMatrix[tour.get(tour.size() - 1)][tour.get(0)];
        return cost;
    }
    public static int calculateDistance(double x1, double y1, double x2, double y2) {
        return (int) ((Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2))) * 10000);
    }

    public static void main(String[] args) throws FileNotFoundException {
        File inputFile = new File("");
        Scanner scanner = new Scanner(inputFile);
        String[] inputconst = scanner.nextLine().split(" ");
        double constant;
        int[][] costMatrix;
        int count = 0;
        if (inputconst.length == 1) {
            constant = Double.parseDouble(inputconst[0]);
            inputconst = scanner.nextLine().split(" ");
            if (inputconst[0].charAt(0) == '(') {
                ArrayList<String> input2 = new ArrayList<>();
                input2.add(inputconst[0] + " " + inputconst[1]);
                while (scanner.hasNext()) {
                    input2.add(scanner.nextLine().replace("\n", ""));
                }
                costMatrix = buildCostMatrix2(input2, constant);
                count = costMatrix.length;
            } else {
                ArrayList<double[]> input2 = new ArrayList<>();
                double[] cordinate1 = new double[2];
                String[] temp;
                cordinate1[0] = Double.parseDouble(inputconst[0]);
                cordinate1[1] = Double.parseDouble(inputconst[1]);
                input2.add(cordinate1);
                while (scanner.hasNext()) {
                    temp = scanner.nextLine().replace("\n", "").split(" ");
                    double[] cordinate = new double[2];
                    cordinate[0] = Double.parseDouble(temp[0]);
                    cordinate[1] = Double.parseDouble(temp[1]);

                    input2.add(cordinate);
                }
                costMatrix = buildCostMatrixFromCoordinates(input2, constant);
                count = costMatrix.length;
            }
        } else {
            List<int[]> edges = new ArrayList<>();
            edges.add(new int[]{Integer.parseInt(inputconst[0]), Integer.parseInt(inputconst[1]), (int) (Double.parseDouble(inputconst[2]) * 10000)});
            count = 1;
            while (scanner.hasNext()) {
                int city1 = scanner.nextInt();
                int city2 = scanner.nextInt();
                int cost = (int) (scanner.nextDouble() * 10000);
                count++;

                edges.add(new int[]{city1, city2, cost});
            }
            int discriminant = 1 + 8 * count;
            int n1 = (int) ((1 + Math.sqrt(discriminant)) / 2);
            int n2 = (int) ((1 - Math.sqrt(discriminant)) / 2);
            count = Math.max(n1, n2);
            costMatrix = buildCostMatrix(edges, count);
        }
        TSPChristofides tsp = new TSPChristofides();
        List<Integer> tour = tsp.christofidesAlgorithm(costMatrix);
        int minCost = tsp.calculateTourCost(tour.subList(0, tour.size() - 1), costMatrix);

        System.out.println("Approximate TSP tour:");
        for (int city : tour) {
            System.out.print((city + 1) + " ");
        }
        System.out.println();
        System.out.println("Minimum cost: " + minCost / 10000.0);
    }
    public ArrayList<Integer> christofidesAlgorithm(int[][] costMatrix) {
        int numberOfCities = costMatrix.length;

        List<int[]> mst = findMST(costMatrix);
        List<Integer> oddVertices = findOddDegreeVertices(mst, numberOfCities);
        List<int[]> matching = findMinimumWeightMatching(oddVertices, costMatrix);
        List<int[]> multigraph = new ArrayList<>(mst);
        multigraph.addAll(matching);

        List<Integer> eulerianCircuit = findEulerianCircuit(multigraph, numberOfCities);
        return convertToHamiltonianCircuit(eulerianCircuit);
    }


}
