import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class TSPConvexHull {
    public static void main(String[] args) throws FileNotFoundException {
        ArrayList<Point2D> coordinates = new ArrayList<>();
        File inputFile = new File("");
        Scanner fileScanner = new Scanner(inputFile);
        String[] firstLine = fileScanner.nextLine().split(" ");
        double scalingFactor = Double.parseDouble(firstLine[0]);
        while (fileScanner.hasNext()) {
            String[] line = fileScanner.nextLine().trim().split(" ");
            coordinates.add(new Point2D.Double(Double.parseDouble(line[0]), Double.parseDouble(line[1])));
        }

        ArrayList<Point2D> optimizedRoute = calculateOptimalRoute(coordinates);
        System.out.println("Optimized Route:");
        for (Point2D point : optimizedRoute) {
            System.out.println(point);
        }
        int totalCost = computeRouteCost(optimizedRoute);
        System.out.println("Total Cost (Scaled): " + totalCost * scalingFactor);
    }
    public static int computeRouteCost(ArrayList<Point2D> route) {
        double totalDistance = 0.0;

        for (int i = 0; i < route.size() - 1; i++) {
            Point2D current = route.get(i);
            Point2D next = route.get(i + 1);
            totalDistance += calculateEuclideanDistance(current, next);
        }
        totalDistance += calculateEuclideanDistance(route.get(route.size() - 1), route.get(0));

        return (int) Math.round(totalDistance);
    }
    public static double calculateEuclideanDistance(Point2D p1, Point2D p2) {
        return Math.sqrt(Math.pow(p2.getX() - p1.getX(), 2) + Math.pow(p2.getY() - p1.getY(), 2));
    }
    public static ArrayList<Point2D> calculateOptimalRoute(ArrayList<Point2D> coordinates) {
        ArrayList<Point2D> hullPoints = findConvexHull(coordinates);
        ArrayList<Point2D> remainingCoordinates = new ArrayList<>(coordinates);
        remainingCoordinates.removeAll(hullPoints);
        ArrayList<Point2D> route = new ArrayList<>(hullPoints);
        for (Point2D point : remainingCoordinates) {
            insertPointGreedily(route, point);
        }

        return route;
    }
    private static ArrayList<Point2D> findConvexHull(ArrayList<Point2D> coordinates) {
        if (coordinates.size() <= 3) {
            return new ArrayList<>(coordinates);
        }
        Point2D anchor = Collections.min(coordinates, Comparator.comparingDouble(Point2D::getY)
                .thenComparingDouble(Point2D::getX));
        ArrayList<Point2D> sortedCoordinates = new ArrayList<>(coordinates);
        sortedCoordinates.sort((p1, p2) -> {
            double angle1 = Math.atan2(p1.getY() - anchor.getY(), p1.getX() - anchor.getX());
            double angle2 = Math.atan2(p2.getY() - anchor.getY(), p2.getX() - anchor.getX());
            return Double.compare(angle1, angle2);
        });
        ArrayList<Point2D> hull = new ArrayList<>();
        for (Point2D point : sortedCoordinates) {
            while (hull.size() >= 2 && computeCrossProduct(hull.get(hull.size() - 2), hull.get(hull.size() - 1), point) <= 0) {
                hull.remove(hull.size() - 1);
            }
            hull.add(point);
        }

        return hull;
    }
    private static double computeCrossProduct(Point2D p1, Point2D p2, Point2D p3) {
        double dx1 = p2.getX() - p1.getX();
        double dy1 = p2.getY() - p1.getY();
        double dx2 = p3.getX() - p2.getX();
        double dy2 = p3.getY() - p2.getY();
        return dx1 * dy2 - dy1 * dx2;
    }
    private static void insertPointGreedily(ArrayList<Point2D> route, Point2D point) {
        int bestInsertionIndex = 0;
        double minimumIncrease = Double.MAX_VALUE;

        for (int i = 0; i < route.size(); i++) {
            Point2D current = route.get(i);
            Point2D next = route.get((i + 1) % route.size());
            double increase = calculateEuclideanDistance(current, point)
                    + calculateEuclideanDistance(point, next)
                    - calculateEuclideanDistance(current, next);

            if (increase < minimumIncrease) {
                minimumIncrease = increase;
                bestInsertionIndex = i + 1;
            }
        }
        route.add(bestInsertionIndex, point);
    }
}
