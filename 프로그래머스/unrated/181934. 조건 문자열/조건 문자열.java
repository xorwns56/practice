class Solution {
    public int solution(String ineq, String eq, int n, int m) {
        if((ineq.equals("<") && n < m) || (ineq.equals(">") && n > m) || (eq.equals("=") && n == m)) return 1;
        return 0;
    }
}