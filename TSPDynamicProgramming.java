import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPDynamicProgramming {

    public static ArrayList<Integer> tspHeldKarp(double[][] costMatrix) {
        int n = costMatrix.length;
        int allStates = 1 << n;
        double[][] dp = new double[allStates][n];

        for (double[] row : dp) {
            Arrays.fill(row, Double.MAX_VALUE / 2);
        }
        dp[1][0] = 0;

        for (int state = 0; state < allStates; state++) {
            for (int currentCity = 0; currentCity < n; currentCity++) {
                if ((state & (1 << currentCity)) == 0) continue;

                for (int nextCity = 0; nextCity < n; nextCity++) {
                    if ((state & (1 << nextCity)) != 0) continue;

                    int nextState = state | (1 << nextCity);
                    dp[nextState][nextCity] = Math.min(
                            dp[nextState][nextCity],
                            dp[state][currentCity] + costMatrix[currentCity][nextCity]
                    );
                }
            }
        }

        double minCost = Double.MAX_VALUE;
        int finalState = allStates - 1;
        for (int lastCity = 1; lastCity < n; lastCity++) {
            minCost = Math.min(minCost, dp[finalState][lastCity] + costMatrix[lastCity][0]);
        }

        ArrayList<Integer> tour = new ArrayList<>();
        int lastIndex = 0;
        int state = finalState;
        tour.add(0);

        for (int i = 1; i < n; i++) {

            int bestIndex = -1;
            double bestDist = Double.MAX_VALUE;
            for (int j = 0; j < n; j++) {
                if (j == 0 || notIn(j, state)) continue;

                double newDist = dp[state][j] + costMatrix[j][lastIndex];
                if (newDist < bestDist) {
                    bestIndex = j;
                    bestDist = newDist;
                }
            }

            tour.add(bestIndex);
            state = state ^ (1 << bestIndex);
            lastIndex = bestIndex;
        }

        tour.add(0);

        System.out.println("Minimum Cost: " + minCost);
        System.out.println("Optimal Tour: " + tour);

        return tour;
    }

    private static boolean notIn(int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    public static double calculateDistance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }

    public static double[][] buildCostMatrix(List<double[]> edges, int numberOfCities) {
        double[][] costMatrix = new double[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Double.MAX_VALUE / 2);
        }

        for (double[] edge : edges) {
            double city1 = edge[0];
            double city2 = edge[1];
            double cost = edge[2];
            costMatrix[(int) city1][(int) city2] = cost;
            costMatrix[(int) city2][(int)city1] = cost;
        }

        return costMatrix;
    }

    public static double[][] buildCostMatrix2(List<String> edgesInput ,double constant) {
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
        double[][] costMatrix = new double[numberOfCities][numberOfCities];
        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Double.MAX_VALUE);
        }

        for (double[] edge : edges) {
            int city1 = (int) edge[0];
            int city2 = (int) edge[1];
            double distance = edge[2];
            costMatrix[city1][city2] = constant * distance;
            costMatrix[city2][city1] = constant * distance;
        }

        return costMatrix;
    }

    public static double[][] buildCostMatrixFromCoordinates(List<double[]> coordinates,double cost) {
        int numberOfCities = coordinates.size();
        double[][] costMatrix = new double[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Double.MAX_VALUE);
        }

        for (int i = 0; i < numberOfCities; i++) {
            for (int j = 0; j < numberOfCities; j++) {
                if (i != j) {
                    double[] point1 = coordinates.get(i);
                    double[] point2 = coordinates.get(j);
                    double distance = calculateDistance(point1[0], point1[1], point2[0], point2[1]);
                    costMatrix[i][j] = distance*cost;
                    costMatrix[j][i] = distance*cost;
                } else {
                    costMatrix[i][j] = 0;
                }
            }
        }
        return costMatrix;
    }

    public static void main(String[] args) throws FileNotFoundException {
        File inputFile = new File("");
        Scanner scanner = new Scanner(inputFile);
        String[] inputconst=scanner.nextLine().split(" ");
        double constant;
        double[][] costMatrix;
        int count=0;
        if (inputconst.length==1){
            constant=Integer.parseInt(inputconst[0]);
            inputconst=scanner.nextLine().split(" ");
            if (inputconst[0].charAt(0)=='('){
                ArrayList<String> input2=new ArrayList<>();
                input2.add(inputconst[0]+" "+inputconst[1]);
                while (scanner.hasNext()){
                    input2.add(scanner.nextLine().replace("\n",""));
                }
                costMatrix=buildCostMatrix2(input2,constant);
                count=costMatrix.length;

            }else {
                ArrayList<double[]> input2=new ArrayList<>();
                double[] cordinate1=new double[2];
                String[] temp;
                cordinate1[0]=Double.parseDouble(inputconst[0]);
                cordinate1[1]=Double.parseDouble(inputconst[1]);
                input2.add(cordinate1);
                while (scanner.hasNext()){
                    temp=scanner.nextLine().replace("\n","").split(" ");
                    double[] cordinate=new double[2];
                    cordinate[0]=Double.parseDouble(temp[0]);
                    cordinate[1]=Double.parseDouble(temp[1]);

                    input2.add(cordinate);
                }
                costMatrix=buildCostMatrixFromCoordinates(input2,constant);
                count=costMatrix.length;
            }
        }
        else {
            List<double[]> edges = new ArrayList<>();
            edges.add(new double []{Integer.parseInt(inputconst[0]), Integer.parseInt(inputconst[1]), Double.parseDouble(inputconst[2])});
            count=1;
            while (scanner.hasNext()) {
                int  city1 = scanner.nextInt();
                int city2 = scanner.nextInt();
                double cost = scanner.nextDouble();
                count++;

                edges.add(new double []{city1, city2, cost});
            }
            int discriminant = 1 + 8 * count;
            int n1 = (int) ((1 + Math.sqrt(discriminant)) / 2);
            int n2 = (int) ((1 - Math.sqrt(discriminant)) / 2);
            count=Math.max(n1,n2);
            costMatrix = buildCostMatrix(edges,count);
        }
        tspHeldKarp(costMatrix);

        scanner.close();
    }
}

