//! Tests of the DIA retention-time calibration sizes.

use super::calibrate::dia_rt_cal_sizes;

#[test]
fn the_rt_calibration_sample_doubles_up_to_the_whole_run() {
    assert_eq!(
        dia_rt_cal_sizes(292_883),
        [5_858, 11_716, 23_432, 46_864, 93_728, 187_456, 292_883]
    );
    assert_eq!(dia_rt_cal_sizes(12_032), [2_000, 4_000, 8_000, 12_032]);
    assert_eq!(dia_rt_cal_sizes(1_500), [1_500]);
    assert_eq!(dia_rt_cal_sizes(4_000), [2_000, 4_000]);
}
