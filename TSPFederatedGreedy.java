import java.io.File;
import java.io.IOException;
import java.util.*;

public class TSPFederatedGreedy {

    static class City {
        int id;
        Map<Integer, Double> distances;

        public City(int id) {
            this.id = id;
            this.distances = new HashMap<>();
        }

        public void addDistance(int toCity, double distance) {
            distances.put(toCity, distance);
        }

        public int getNearestUnvisitedCity(Set<Integer> visited) {
            double minDistance = Double.MAX_VALUE;
            int nearestCity = -1;
            for (Map.Entry<Integer, Double> entry : distances.entrySet()) {
                int city = entry.getKey();
                double distance = entry.getValue();
                if (!visited.contains(city) && distance < minDistance) {
                    minDistance = distance;
                    nearestCity = city;
                }
            }
            return nearestCity;
        }
    }

    static class FederatedServer {
        List<City> cities;
        List<Integer> globalRoute;

        public FederatedServer(List<City> cities) {
            this.cities = cities;
            this.globalRoute = new ArrayList<>();
        }

        public List<Integer> solveTSP(int startCity) {
            Set<Integer> visited = new HashSet<>();
            int currentCity = startCity;
            globalRoute.add(currentCity);
            visited.add(currentCity);

            while (visited.size() < cities.size()) {
                int nextCity = -1;
                double minDistance = Double.MAX_VALUE;

                for (City city : cities) {
                    if (city.id == currentCity) {
                        int localNextCity = city.getNearestUnvisitedCity(visited);
                        if (localNextCity != -1 && city.distances.get(localNextCity) < minDistance) {
                            minDistance = city.distances.get(localNextCity);
                            nextCity = localNextCity;
                        }
                    }
                }

                if (nextCity == -1) {
                    throw new RuntimeException("No valid route found");
                }

                globalRoute.add(nextCity);
                visited.add(nextCity);
                currentCity = nextCity;
            }

            globalRoute.add(startCity);
            return globalRoute;
        }

        public double calculateRouteCost(List<Integer> route) {
            double totalCost = 0.0;
            for (int i = 0; i < route.size() - 1; i++) {
                int fromCity = route.get(i);
                int toCity = route.get(i + 1);
                totalCost += cities.get(fromCity).distances.get(toCity);
            }
            return totalCost;
        }
    }

    public static void main(String[] args) throws IOException {
        File inputFile = new File("");
        Scanner scanner = new Scanner(inputFile);
        String[] initialInput = scanner.nextLine().split(" ");
        double scaleFactor;
        int[][] distanceMatrix;
        int cityCount = 0;

        if (initialInput.length == 1) {
            scaleFactor = Double.parseDouble(initialInput[0]);
            initialInput = scanner.nextLine().split(" ");
            if (initialInput[0].charAt(0) == '(') {
                ArrayList<String> edgeList = new ArrayList<>();
                edgeList.add(initialInput[0] + " " + initialInput[1]);
                while (scanner.hasNext()) {
                    edgeList.add(scanner.nextLine().replace("\n", ""));
                }
                distanceMatrix = createDistanceMatrixFromEdges(edgeList, scaleFactor);
                cityCount = distanceMatrix.length;
            } else {
                ArrayList<double[]> coordinatesList = new ArrayList<>();
                double[] coord1 = new double[2];
                String[] temp;
                coord1[0] = Double.parseDouble(initialInput[0]);
                coord1[1] = Double.parseDouble(initialInput[1]);
                coordinatesList.add(coord1);
                while (scanner.hasNext()) {
                    temp = scanner.nextLine().replace("\n", "").split(" ");
                    double[] coord = new double[2];
                    coord[0] = Double.parseDouble(temp[0]);
                    coord[1] = Double.parseDouble(temp[1]);

                    coordinatesList.add(coord);
                }
                distanceMatrix = createDistanceMatrixFromCoordinates(coordinatesList, scaleFactor);
                cityCount = distanceMatrix.length;
            }
        } else {
            List<int[]> edgeList = new ArrayList<>();
            edgeList.add(new int[]{Integer.parseInt(initialInput[0]), Integer.parseInt(initialInput[1]), (int) (Double.parseDouble(initialInput[2]) * 10000)});
            cityCount = 1;
            while (scanner.hasNext()) {
                int city1 = scanner.nextInt();
                int city2 = scanner.nextInt();
                int cost = (int) (scanner.nextDouble() * 10000);
                cityCount++;

                edgeList.add(new int[]{city1, city2, cost});
            }
            int discriminant = 1 + 8 * cityCount;
            int n1 = (int) ((1 + Math.sqrt(discriminant)) / 2);
            int n2 = (int) ((1 - Math.sqrt(discriminant)) / 2);
            cityCount = Math.max(n1, n2);
            distanceMatrix = createDistanceMatrixFromEdges(edgeList, cityCount);
        }

        FederatedServer server = new FederatedServer(createCities(distanceMatrix));
        List<Integer> tour = server.solveTSP(0);
        double cost = server.calculateRouteCost(tour);

        System.out.println("Optimal Tour: " + tour);
        System.out.println("Minimum Cost: " + cost / 10000.0);
    }



    public static List<City> createCities(int[][] distanceMatrix) {
        List<City> cities = new ArrayList<>();
        for (int i = 0; i < distanceMatrix.length; i++) {
            City city = new City(i);
            for (int j = 0; j < distanceMatrix[i].length; j++) {
                if (i != j && distanceMatrix[i][j] != Integer.MAX_VALUE) {
                    city.addDistance(j, distanceMatrix[i][j]);
                }
            }
            cities.add(city);
        }
        return cities;
    }

    public static int[][] createDistanceMatrixFromEdges(List<int[]> edges, int numberOfCities) {
        int[][] matrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(matrix[i], Integer.MAX_VALUE / 2);
        }

        for (int[] edge : edges) {
            int city1 = edge[0];
            int city2 = edge[1];
            int cost = edge[2];
            matrix[city1][city2] = cost;
            matrix[city2][city1] = cost;
        }

        return matrix;
    }

    public static int[][] createDistanceMatrixFromCoordinates(List<double[]> coordinates, double scaleFactor) {
        int numCities = coordinates.size();
        int[][] matrix = new int[numCities][numCities];

        for (int i = 0; i < numCities; i++) {
            Arrays.fill(matrix[i], Integer.MAX_VALUE);
        }

        for (int i = 0; i < numCities; i++) {
            for (int j = 0; j < numCities; j++) {
                if (i != j) {
                    double[] city1 = coordinates.get(i);
                    double[] city2 = coordinates.get(j);
                    double distance = calculateEuclideanDistance(city1[0], city1[1], city2[0], city2[1]);
                    matrix[i][j] = (int) (distance * scaleFactor);
                    matrix[j][i] = (int) (distance * scaleFactor);
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
        return matrix;
    }

    public static int[][] createDistanceMatrixFromEdges(List<String> edgeInput, double scaleFactor) {
        Map<Integer, double[]> cityCoordinates = new HashMap<>();
        Map<String, Integer> coordinateToIndex = new HashMap<>();

        int cityIndex = 0;
        List<double[]> edges = new ArrayList<>();
        for (String edge : edgeInput) {
            String[] points = edge.split(" ");
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

            double distance = calculateEuclideanDistance(point1[0], point1[1], point2[0], point2[1]);

            edges.add(new double[]{city1, city2, distance});
        }

        int numCities = cityIndex;
        int[][] matrix = new int[numCities][numCities];
        for (int i = 0; i < numCities; i++) {
            Arrays.fill(matrix[i], Integer.MAX_VALUE);
        }

        for (double[] edge : edges) {
            int city1 = (int) edge[0];
            int city2 = (int) edge[1];
            double distance = edge[2];
            matrix[city1][city2] = (int) (scaleFactor * distance);
            matrix[city2][city1] = (int) (scaleFactor * distance);
        }

        return matrix;
    }

    public static int calculateEuclideanDistance(double x1, double y1, double x2, double y2) {
        return (int) ((Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2))) * 10000);
    }


}