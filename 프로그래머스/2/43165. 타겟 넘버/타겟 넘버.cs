using System;
public class Solution {
    public int solution(int[] numbers, int target) {
        return Cal(numbers, -1, 0, target);
    }
    int Cal(int[] numbers, int i, int sum, int target){
        if(i == numbers.Length - 1) return sum == target ? 1 : 0;
        return Cal(numbers, i + 1, sum + numbers[i + 1], target)
            + Cal(numbers, i + 1, sum - numbers[i + 1], target);
    }
}