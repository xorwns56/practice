class Solution {
    public int solution(int[] dot) {
        return (dot[0] * dot[1] < 0 ? 2 : 1) + (dot[1] < 0 ? 2 : 0);
    }
}