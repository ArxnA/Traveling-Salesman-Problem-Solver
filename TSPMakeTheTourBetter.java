import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPMakeTheTourBetter {

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
    public static int[][] buildCostMatrix2(List<String> edgesInput ,double constant) {
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
            costMatrix[city1][city2] = (int)(constant * distance);
            costMatrix[city2][city1] = (int)( constant * distance);
        }

        return costMatrix;
    }
    public static int[][] buildCostMatrixFromCoordinates(List<double[]> coordinates,double cost) {
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
                    costMatrix[i][j] = (int)(distance*cost);
                    costMatrix[j][i] = (int)(distance*cost);
                } else {
                    costMatrix[i][j] = 0;
                }
            }
        }
        return costMatrix;
    }

    public static int calculateDistance(double x1, double y1, double x2, double y2) {
        return (int)((Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2)))*10000);
    }
    public int calculateTourCost(List<Integer> tour, int[][] costMatrix) {
        int cost = 0;
        for (int i = 0; i < tour.size() - 1; i++) {
            cost += costMatrix[tour.get(i)][tour.get(i + 1)];
        }
        cost += costMatrix[tour.get(tour.size() - 1)][tour.get(0)];
        return cost;
    }

    private List<Integer> twoOptSwap(List<Integer> tour, int i, int k) {
        List<Integer> newTour = new ArrayList<>(tour.subList(0, i));
        List<Integer> reversedSegment = new ArrayList<>(tour.subList(i, k + 1));
        Collections.reverse(reversedSegment);
        newTour.addAll(reversedSegment);
        newTour.addAll(tour.subList(k + 1, tour.size()));
        return newTour;
    }

    public static void main(String[] args) throws FileNotFoundException {
        File inputFile = new File("");
        Scanner scanner = new Scanner(inputFile);
        String[] inputconst=scanner.nextLine().split(" ");
        double constant;
        int[][] costMatrix;
        int count=0;
        if (inputconst.length==1){
            constant=Double.parseDouble(inputconst[0]);
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
        } else {
            List<int[]> edges = new ArrayList<>();
            edges.add(new int []{Integer.parseInt(inputconst[0]), Integer.parseInt(inputconst[1]),(int)(Double.parseDouble(inputconst[2])*10000)});
            count=1;
            while (scanner.hasNext()) {
                int  city1 = scanner.nextInt();
                int city2 = scanner.nextInt();
                int cost = (int)(scanner.nextDouble()*10000);
                count++;
                edges.add(new int []{city1, city2, cost});
            }
            int discriminant = 1 + 8 * count;
            int n1 = (int) ((1 + Math.sqrt(discriminant)) / 2);
            int n2 = (int) ((1 - Math.sqrt(discriminant)) / 2);
            count=Math.max(n1,n2);
            costMatrix = buildCostMatrix(edges,count);
        }
        TSPMakeTheTourBetter tsp = new TSPMakeTheTourBetter();
        List<Integer> tour = tsp.MakeTheTourBetterAlgorithm(costMatrix);
        int minCost = tsp.calculateTourCost(tour, costMatrix);

        System.out.println("Approximate TSP tour:");
        for (int city : tour) {
            System.out.print((city + 1) + " ");
        }
        System.out.println();
        System.out.println("Minimum cost: " + minCost / 10000.0);
    }
    public List<Integer> MakeTheTourBetterAlgorithm(int[][] costMatrix) {
        int numberOfCities = costMatrix.length;

        List<Integer> currentTour = new ArrayList<>();
        for (int i = 0; i < numberOfCities; i++) {
            currentTour.add(i);
        }
        //Collections.shuffle(currentTour);

        List<Integer> bestTour = new ArrayList<>(currentTour);

        boolean improved = true;

        while (improved) {
            improved = false;

            for (int i = 1; i < numberOfCities - 1; i++) {
                for (int k = i + 1; k < numberOfCities; k++) {
                    if (i == 0 && k == numberOfCities - 1) continue;
                    int bestDistance = costMatrix[bestTour.get(i - 1)][bestTour.get(i)] + costMatrix[bestTour.get(k - 1)][bestTour.get(k)];
                    int newDistance = costMatrix[bestTour.get(i - 1)][bestTour.get(k - 1)] + costMatrix[bestTour.get(i)][bestTour.get(k)];

                    if (newDistance < bestDistance) {
                        Collections.reverse(bestTour.subList(i, k));
                        improved = true;
                        break;
                    }
                }
                if (improved) break;
            }
        }

        return bestTour;
    }
}
