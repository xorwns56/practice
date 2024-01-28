import java.util.*;
class Solution {
    public int[] solution(int[] numlist, int n) {
        return Arrays.stream(numlist).boxed().sorted((a, b)->Math.abs(a - n) == Math.abs(b - n) ? b - a : Math.abs(a - n) - Math.abs(b - n)).mapToInt(Integer::valueOf).toArray();
    }
}