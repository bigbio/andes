use super::RefTracker;

#[test]
fn competitor_excludes_the_rows_own_peptidoform() {
    let mut t = RefTracker::default();
    t.fold(30, 1);
    t.fold(25, 2);
    t.fold(28, 1);
    assert_eq!(t.best_excluding(1), Some(25));
    assert_eq!(t.best_excluding(2), Some(30));
    assert_eq!(t.best_excluding(3), Some(30));
    t.fold(40, 3);
    assert_eq!(t.best_excluding(3), Some(30));
    assert_eq!(RefTracker::default().best_excluding(1), None);
    let mut one = RefTracker::default();
    one.fold(30, 1);
    assert_eq!(one.best_excluding(1), None);
    assert_eq!(one.best_excluding(2), Some(30));
}
