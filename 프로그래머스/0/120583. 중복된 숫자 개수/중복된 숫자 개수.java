import java.util.*;
class Solution {
    public int solution(int[] array, int n) {
        return (int)Arrays.stream(array).filter(x->x==n).count();
    }
}