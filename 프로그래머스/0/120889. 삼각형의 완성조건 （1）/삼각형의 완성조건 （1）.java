class Solution {
    public int solution(int[] sides) {
        int max = 0;
        if(sides[max] < sides[1]) max = 1;
        if(sides[max] < sides[2]) max = 2;
        return sides[max] < sides[(max + 1) % 3] + sides[(max + 2) % 3] ? 1 : 2;
    }
}