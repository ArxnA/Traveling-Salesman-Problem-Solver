import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPBruteForce {

    public static double calculateTourCost(double[][] costMatrix, int[] tour) {
        double totalCost = 0;
        for (int i = 0; i < tour.length - 1; i++) {
            totalCost += costMatrix[tour[i]][tour[i + 1]];
        }
        totalCost += costMatrix[tour[tour.length - 1]][tour[0]];
        return totalCost;
    }

    public static ArrayList<Integer> findMinimumCost(double[][] costMatrix, int numberOfCities) {
        int[] cities = new int[numberOfCities];
        for (int i = 0; i < numberOfCities; i++) {
            cities[i] = i;
        }

        double minCost = Double.MAX_VALUE;
        int[] bestTour = null;

        do {
            double currentCost = calculateTourCost(costMatrix, cities);
            if (currentCost < minCost) {
                minCost = currentCost;
                bestTour = cities.clone();
            }
        } while (nextPermutation(cities));

        ArrayList<Integer> tour1=new ArrayList<>();
        System.out.println("Minimum Cost: " + minCost);
        System.out.print("Best Tour: ");
        if (bestTour != null) {
            for (int i = 0; i <bestTour.length ; i++) {
                System.out.print(bestTour[i] + " ");
                tour1.add(bestTour[i]);
            }
            System.out.println(bestTour[0]);
        }

        return tour1;
    }

    public static boolean nextPermutation(int[] array) {
        int i = array.length - 2;
        while (i >= 0 && array[i] >= array[i + 1]) {
            i--;
        }
        if (i < 0) return false;

        int j = array.length - 1;
        while (array[j] <= array[i]) {
            j--;
        }

        swap(array, i, j);

        int left = i + 1, right = array.length - 1;
        while (left < right) {
            swap(array, left++, right--);
        }

        return true;
    }

    public static void swap(int[] array, int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    public static double[][] buildCostMatrix(List<double[]> edges, int numberOfCities) {
        double[][] costMatrix = new double[numberOfCities][numberOfCities];
        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Double.MAX_VALUE);
        }
        for (double[] edge : edges) {
            double city1 = edge[0] ;
            double city2 = edge[1] ;
            double cost = edge[2];
            costMatrix[(int)city1][(int) city2] = cost;
            costMatrix[(int) city2][(int) city1] = cost;
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

    public static double calculateDistance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
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
        } else {
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

        findMinimumCost(costMatrix, count);
        scanner.close();
    }
}
