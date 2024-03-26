import java.util.*;
class Solution {
    public String solution(int[] numbers) {
        String[] numbers_str = new String[numbers.length];
        for(int i = 0; i < numbers.length; i++) numbers_str[i] = Integer.toString(numbers[i]);
        Arrays.sort(numbers_str, (s1, s2)->{ return (s2 + s1).compareTo(s1 + s2); });
        if(numbers_str[0].equals("0")) return "0";
        StringBuilder stringBuilder = new StringBuilder();
        for(int i = 0; i < numbers_str.length; i++) stringBuilder.append(numbers_str[i]);
        return stringBuilder.toString();
    }
}