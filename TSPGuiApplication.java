import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;
import javafx.stage.Stage;
import java.util.*;
import java.util.stream.IntStream;

public class TSPGuiApplication extends Application {

    private List<int[]> cities = new ArrayList<>();
    private List<Integer> tour = new ArrayList<>();
    private final Canvas canvas = new Canvas(1200, 600);
    private Label tourCostLabel = new Label("Tour Cost: N/A");
    private int flag;

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) {
        StackPane root = new StackPane();
        root.getChildren().add(canvas);

        tourCostLabel.setTranslateY(280);
        tourCostLabel.setTranslateX(-550);
        root.getChildren().add(tourCostLabel);

        Button dynamicProgrammingButton = new Button("Run Dynamic Programming");
        Button christofidesButton = new Button("Run Christofides");
        Button bruteForceButton = new Button("Run Brute Force");
        Button shortestNextButton = new Button("Run Shortest Next");
        Button optimizeShortestNextButton = new Button("Run Optimize Shortest Next");
        Button disjointCycleButton = new Button("Run Disjoint Cycle");
        Button makeTheTourBetterButton = new Button("Run Make The Tour Better");
        Button bestPathFromPathsButton = new Button("Run Best Path From Paths");
        Button shortestEdgeFirstButton = new Button("Run Shortest Edge First");
        Button optimizeShortestEdgeFirstButton = new Button("Run Optimize Shortest Edge First");
        Button federatedGreedyButton = new Button("Run Federated Greedy");
        Button clearPointsButton = new Button("Clear Points");


        clearPointsButton.setTranslateY(280);
        clearPointsButton.setTranslateX(550);

        dynamicProgrammingButton.setTranslateY(-270);
        dynamicProgrammingButton.setTranslateX(-490);

        christofidesButton.setTranslateY(-230);
        christofidesButton.setTranslateX(-490);

        bruteForceButton.setTranslateY(-190);
        bruteForceButton.setTranslateX(-490);

        shortestNextButton.setTranslateY(-150);
        shortestNextButton.setTranslateX(-490);

        optimizeShortestNextButton.setTranslateY(-110);
        optimizeShortestNextButton.setTranslateX(-490);

        disjointCycleButton.setTranslateY(-70);
        disjointCycleButton.setTranslateX(-490);

        makeTheTourBetterButton.setTranslateY(-30);
        makeTheTourBetterButton.setTranslateX(-490);

        bestPathFromPathsButton.setTranslateY(10);
        bestPathFromPathsButton.setTranslateX(-490);

        shortestEdgeFirstButton.setTranslateY(50);
        shortestEdgeFirstButton.setTranslateX(-490);

        optimizeShortestEdgeFirstButton.setTranslateY(90);
        optimizeShortestEdgeFirstButton.setTranslateX(-490);

        federatedGreedyButton.setTranslateY(130);
        federatedGreedyButton.setTranslateX(-490);


        root.getChildren().addAll(christofidesButton, dynamicProgrammingButton, bruteForceButton, shortestNextButton, disjointCycleButton, makeTheTourBetterButton, bestPathFromPathsButton, optimizeShortestNextButton, shortestEdgeFirstButton, optimizeShortestEdgeFirstButton, federatedGreedyButton, clearPointsButton);

        christofidesButton.setOnAction(event -> runChristofidesTSP());
        dynamicProgrammingButton.setOnAction(event -> runDynamicProgrammingTSP());
        bruteForceButton.setOnAction(event -> runBruteForceTSP());
        shortestNextButton.setOnAction(event -> runShortestNextTSP());
        disjointCycleButton.setOnAction(event -> runDisjointCycleTSP());
        makeTheTourBetterButton.setOnAction(event -> runMakeTheTourBetterTSP());
        bestPathFromPathsButton.setOnAction(event -> runBestPathFromPathsTSP());
        optimizeShortestNextButton.setOnAction(event -> runOptimizeShortestNextTSP());
        shortestEdgeFirstButton.setOnAction(event -> runShortestEdgeFirstTSP());
        optimizeShortestEdgeFirstButton.setOnAction(event -> runShortestEdgeFirstOptimizeTSP());
        federatedGreedyButton.setOnAction(event -> runFederatedGreedyTSP());

        clearPointsButton.setOnAction(event -> {
            cities.clear();
            tour.clear();
            drawCities();
        });

        canvas.setOnMouseClicked(event -> addCity(event.getX(), event.getY()));

        primaryStage.setTitle("TSP Algorithms");
        primaryStage.setScene(new Scene(root, 1200, 600));
        primaryStage.show();
    }

    private void addCity(double x, double y) {
        cities.add(new int[]{(int) x, (int) y});
        drawCities();
    }

    private void drawCities() {
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());

        gc.setFill(Color.BLUE);
        for (int[] city : cities) {
            gc.fillOval(city[0] - 5, city[1] - 5, 10, 10);
        }

        if (!tour.isEmpty()) {
            gc.setStroke(Color.RED);
            for (int i = 0; i < tour.size() - 1; i++) {
                int start = tour.get(i);
                int end = tour.get(i + 1);
                gc.strokeLine(cities.get(start)[0], cities.get(start)[1], cities.get(end)[0], cities.get(end)[1]);
            }
            int start = tour.get(tour.size() - 1);
            int end = tour.get(0);
            gc.strokeLine(cities.get(start)[0], cities.get(start)[1], cities.get(end)[0], cities.get(end)[1]);
        }

        updateTourCost();
    }

    private void updateTourCost() {
        if (tour.isEmpty() || cities.size() < 2) {
            tourCostLabel.setText("Tour Cost: N/A");
            return;
        }

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        int totalCost = 0;
        if (flag==1){
            totalCost-=1073741823;
        }

        for (int i = 0; i < tour.size() - 1; i++) {
            totalCost += costMatrix[tour.get(i)][tour.get(i + 1)];
        }
        totalCost += costMatrix[tour.get(tour.size() - 1)][tour.get(0)];
        System.out.println(tour);
        System.out.println(costMatrix[tour.get(tour.size() - 1)][tour.get(0)]);

        tourCostLabel.setText("Tour Cost: " + totalCost);
    }


    private void runFederatedGreedyTSP() {
        tour.clear();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        TSPFederatedGreedy.FederatedServer server = new TSPFederatedGreedy.FederatedServer(TSPFederatedGreedy.createCities(costMatrix));
        tour = server.solveTSP(0);
        flag=1;

        drawCities();
    }

    private void runChristofidesTSP() {
        tour.clear();
        TSPChristofides tsp = new TSPChristofides();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.christofidesAlgorithm(costMatrix);
        flag=1;

        drawCities();
    }
    private void runShortestNextTSP() {
        tour.clear();
        TSPOptimizeShortestNext tsp = new TSPOptimizeShortestNext();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.nearestNeighborAlgorithm(costMatrix);
        flag=1;

        drawCities();
    }
    private void runOptimizeShortestNextTSP() {
        tour.clear();
        TSPOptimizeShortestNext tsp = new TSPOptimizeShortestNext();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.nearestNeighborAlgorithm(costMatrix);
        tour = tsp.twoOptOptimization(tour, costMatrix);
        flag=1;

        drawCities();
    }
    private void runShortestEdgeFirstTSP() {
        tour.clear();
        TSPShortestEdgeFirst tsp = new TSPShortestEdgeFirst();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.greedyTSP(costMatrix);
        flag=0;

        drawCities();
    }
    private void runShortestEdgeFirstOptimizeTSP() {
        tour.clear();
        TSPShortestEdgeFirst tsp = new TSPShortestEdgeFirst();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.greedyTSP(costMatrix);
        tour = tsp.twoOpt(tour, costMatrix);
        flag=0;

        drawCities();
    }
    private void runDisjointCycleTSP() {
        tour.clear();
        TSPDisjointCycle tsp = new TSPDisjointCycle();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.christofidesWithDisjointCycle(costMatrix);
        flag=1;

        drawCities();
    }
    private void runMakeTheTourBetterTSP() {
        tour.clear();
        TSPMakeTheTourBetter tsp = new TSPMakeTheTourBetter();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.MakeTheTourBetterAlgorithm(costMatrix);
        flag=0;

        drawCities();
    }
    private void runBestPathFromPathsTSP() {
        tour.clear();
        TSPBestPathFromPaths tsp = new TSPBestPathFromPaths();

        int[][] costMatrix = buildCostMatrixInt(cities, cities.size());
        tour = tsp.FindTheBestPath(costMatrix);
        flag=0;

        drawCities();
    }

    private void runDynamicProgrammingTSP() {
        tour.clear();
        try {
            double[][] costMatrix = buildCostMatrix(cities, cities.size());
            TSPDynamicProgramming dpTSP = new TSPDynamicProgramming();
            tour = dpTSP.tspHeldKarp(costMatrix);
            flag=1;
            drawCities();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void runBruteForceTSP() {
        tour.clear();
        try {
            double[][] costMatrix = buildCostMatrix(cities, cities.size());
            TSPBruteForce bruteForceTSP = new TSPBruteForce();
            tour = bruteForceTSP.findMinimumCost(costMatrix, cities.size());
            flag=0;
            drawCities();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public static double[][] buildCostMatrix(List<int[]> cities, int numberOfCities) {
        double[][] costMatrix = new double[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Double.MAX_VALUE);
        }

        for (int i = 0; i < cities.size(); i++) {
            for (int j = i + 1; j < cities.size(); j++) {
                double cost = calculateDistance(cities.get(i), cities.get(j));
                costMatrix[i][j] = cost;
                costMatrix[j][i] = cost;
            }
        }

        return costMatrix;
    }

    public static int[][] buildCostMatrixInt(List<int[]> cities, int numberOfCities) {
        System.out.println("Number of cities: " + cities.size());
        int[][] costMatrix = new int[numberOfCities][numberOfCities];

        for (int i = 0; i < numberOfCities; i++) {
            Arrays.fill(costMatrix[i], Integer.MAX_VALUE / 2);
        }

        for (int i = 0; i < cities.size(); i++) {
            for (int j = i + 1; j < cities.size(); j++) {
                int[] city1 = cities.get(i);
                int[] city2 = cities.get(j);
                int cost = calculateDistanceInt(city1, city2);
                costMatrix[i][j] = cost;
                costMatrix[j][i] = cost;
            }
        }

        System.out.println("Cost Matrix:");
        for (int[] row : costMatrix) {
            System.out.println(Arrays.toString(row));
        }

        return costMatrix;
    }

    private static int calculateDistanceInt(int[] city1, int[] city2) {
        return (int) Math.sqrt(Math.pow(city1[0] - city2[0], 2) + Math.pow(city1[1] - city2[1], 2));
    }

    private static double calculateDistance(int[] city1, int[] city2) {
        return Math.sqrt(Math.pow(city1[0] - city2[0], 2) + Math.pow(city1[1] - city2[1], 2));
    }
}
