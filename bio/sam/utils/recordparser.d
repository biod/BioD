module bio.sam.utils.recordparser;

#line 1 "sam_alignment.rl"
/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    BioD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    BioD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#line 24 "sam_alignment.d"
static byte[] _sam_alignment_actions = [
	0, 1, 0, 1, 2, 1, 4, 1, 
	6, 1, 7, 1, 8, 1, 9, 1, 
	10, 1, 11, 1, 12, 1, 15, 1, 
	16, 1, 17, 1, 18, 1, 19, 1, 
	21, 1, 22, 1, 24, 1, 25, 1, 
	27, 1, 31, 1, 35, 1, 36, 1, 
	37, 2, 1, 2, 2, 3, 20, 2, 
	3, 32, 2, 5, 33, 2, 6, 7, 
	2, 13, 14, 2, 23, 24, 2, 29, 
	37, 2, 30, 37, 3, 3, 26, 37, 
	3, 5, 28, 37, 3, 15, 1, 2, 
	4, 3, 32, 34, 37, 4, 5, 33, 
	34, 37
];

static short[] _sam_alignment_key_offsets = [
	0, 5, 6, 9, 10, 17, 18, 21, 
	22, 25, 26, 29, 30, 36, 37, 40, 
	41, 46, 47, 54, 55, 62, 65, 68, 
	71, 74, 77, 80, 83, 86, 89, 92, 
	95, 98, 101, 104, 107, 110, 113, 116, 
	117, 120, 123, 126, 129, 132, 135, 138, 
	141, 144, 147, 150, 153, 156, 159, 162, 
	165, 168, 169, 172, 173, 185, 197, 209, 
	221, 233, 245, 257, 269, 281, 293, 305, 
	317, 329, 341, 353, 365, 377, 387, 390, 
	393, 396, 399, 402, 405, 408, 411, 414, 
	417, 420, 423, 426, 429, 432, 435, 438, 
	441, 442, 445, 448, 451, 454, 457, 460, 
	463, 466, 469, 472, 475, 478, 481, 484, 
	487, 490, 493, 494, 497, 500, 503, 506, 
	509, 512, 515, 518, 521, 524, 527, 530, 
	533, 536, 539, 542, 545, 548, 549, 554, 
	557, 558, 563, 570, 572, 579, 581, 584, 
	585, 587, 595, 597, 602, 605, 609, 613, 
	617, 621, 625, 629, 633, 637, 641, 645, 
	649, 653, 657, 661, 665, 669, 673, 675, 
	677, 685, 690, 693, 699, 704, 707, 711, 
	718, 720, 722, 724, 726, 728, 730, 737, 
	744, 746, 749, 752, 754, 762, 767, 770, 
	775, 780, 783, 786, 792, 794, 796, 797, 
	799, 801, 803, 808, 811, 814, 817, 820, 
	823, 826, 829, 832, 835, 838, 841, 844, 
	847, 850, 853, 856, 859, 862, 863
];

static char[] _sam_alignment_trans_keys = [
	9u, 33u, 63u, 65u, 126u, 9u, 9u, 48u, 
	57u, 9u, 9u, 33u, 41u, 43u, 60u, 62u, 
	126u, 9u, 9u, 48u, 57u, 9u, 9u, 48u, 
	57u, 9u, 9u, 48u, 57u, 9u, 9u, 61u, 
	33u, 41u, 43u, 126u, 9u, 9u, 48u, 57u, 
	9u, 9u, 43u, 45u, 48u, 57u, 9u, 9u, 
	46u, 61u, 65u, 90u, 97u, 122u, 9u, 9u, 
	46u, 61u, 65u, 90u, 97u, 122u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 9u, 33u, 126u, 9u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 48u, 
	57u, 72u, 73u, 77u, 78u, 9u, 61u, 68u, 
	80u, 83u, 88u, 48u, 57u, 72u, 73u, 77u, 
	78u, 9u, 61u, 68u, 80u, 83u, 88u, 72u, 
	73u, 77u, 78u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 9u, 33u, 
	126u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 9u, 33u, 63u, 
	65u, 126u, 9u, 33u, 126u, 9u, 9u, 65u, 
	90u, 97u, 122u, 9u, 48u, 57u, 65u, 90u, 
	97u, 122u, 9u, 58u, 9u, 65u, 66u, 72u, 
	90u, 102u, 105u, 9u, 58u, 9u, 33u, 126u, 
	9u, 9u, 58u, 9u, 67u, 73u, 83u, 99u, 
	102u, 105u, 115u, 9u, 44u, 9u, 43u, 45u, 
	48u, 57u, 9u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 48u, 57u, 9u, 44u, 48u, 
	57u, 9u, 44u, 9u, 44u, 9u, 43u, 45u, 
	46u, 105u, 110u, 48u, 57u, 9u, 46u, 105u, 
	48u, 57u, 9u, 48u, 57u, 9u, 44u, 69u, 
	101u, 48u, 57u, 9u, 43u, 45u, 48u, 57u, 
	9u, 48u, 57u, 9u, 44u, 48u, 57u, 9u, 
	44u, 46u, 69u, 101u, 48u, 57u, 9u, 110u, 
	9u, 102u, 9u, 44u, 9u, 97u, 9u, 110u, 
	9u, 58u, 9u, 48u, 57u, 65u, 70u, 97u, 
	102u, 9u, 48u, 57u, 65u, 70u, 97u, 102u, 
	9u, 58u, 9u, 32u, 126u, 9u, 32u, 126u, 
	9u, 58u, 9u, 43u, 45u, 46u, 105u, 110u, 
	48u, 57u, 9u, 46u, 105u, 48u, 57u, 9u, 
	48u, 57u, 9u, 69u, 101u, 48u, 57u, 9u, 
	43u, 45u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 46u, 69u, 101u, 48u, 57u, 
	9u, 110u, 9u, 102u, 9u, 9u, 97u, 9u, 
	110u, 9u, 58u, 9u, 43u, 45u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 48u, 
	57u, 9u, 48u, 57u, 9u, 48u, 57u, 9u, 
	48u, 57u, 9u, 48u, 57u, 9u, 48u, 57u, 
	9u, 48u, 57u, 9u, 48u, 57u, 9u, 9u, 
	33u, 126u, 0
];

static byte[] _sam_alignment_single_lengths = [
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 2, 1, 1, 1, 
	3, 1, 3, 1, 3, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 6, 6, 6, 6, 
	6, 6, 6, 6, 6, 6, 6, 6, 
	6, 6, 6, 6, 6, 6, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 2, 7, 2, 1, 1, 
	2, 8, 2, 3, 1, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	6, 3, 1, 4, 3, 1, 2, 5, 
	2, 2, 2, 2, 2, 2, 1, 1, 
	2, 1, 1, 2, 6, 3, 1, 3, 
	3, 1, 1, 4, 2, 2, 1, 2, 
	2, 2, 3, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1
];

static byte[] _sam_alignment_range_lengths = [
	2, 0, 1, 0, 3, 0, 1, 0, 
	1, 0, 1, 0, 2, 0, 1, 0, 
	1, 0, 2, 0, 2, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 0, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 0, 1, 0, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 2, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	0, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 0, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 0, 2, 1, 
	0, 2, 3, 0, 0, 0, 1, 0, 
	0, 0, 0, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 0, 0, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	0, 0, 0, 0, 0, 0, 3, 3, 
	0, 1, 1, 0, 1, 1, 1, 1, 
	1, 1, 1, 1, 0, 0, 0, 0, 
	0, 0, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 0, 1
];

static short[] _sam_alignment_index_offsets = [
	0, 4, 6, 9, 11, 16, 18, 21, 
	23, 26, 28, 31, 33, 38, 40, 43, 
	45, 50, 52, 58, 60, 66, 69, 72, 
	75, 78, 81, 84, 87, 90, 93, 96, 
	99, 102, 105, 108, 111, 114, 117, 120, 
	122, 125, 128, 131, 134, 137, 140, 143, 
	146, 149, 152, 155, 158, 161, 164, 167, 
	170, 173, 175, 178, 180, 190, 200, 210, 
	220, 230, 240, 250, 260, 270, 280, 290, 
	300, 310, 320, 330, 340, 350, 359, 362, 
	365, 368, 371, 374, 377, 380, 383, 386, 
	389, 392, 395, 398, 401, 404, 407, 410, 
	413, 415, 418, 421, 424, 427, 430, 433, 
	436, 439, 442, 445, 448, 451, 454, 457, 
	460, 463, 466, 468, 471, 474, 477, 480, 
	483, 486, 489, 492, 495, 498, 501, 504, 
	507, 510, 513, 516, 519, 522, 524, 528, 
	531, 533, 537, 542, 545, 553, 556, 559, 
	561, 564, 573, 576, 581, 584, 588, 592, 
	596, 600, 604, 608, 612, 616, 620, 624, 
	628, 632, 636, 640, 644, 648, 652, 655, 
	658, 666, 671, 674, 680, 685, 688, 692, 
	699, 702, 705, 708, 711, 714, 717, 722, 
	727, 730, 733, 736, 739, 747, 752, 755, 
	760, 765, 768, 771, 777, 780, 783, 785, 
	788, 791, 794, 799, 802, 805, 808, 811, 
	814, 817, 820, 823, 826, 829, 832, 835, 
	838, 841, 844, 847, 850, 853, 855
];

static ubyte[] _sam_alignment_trans_targs = [
	2, 134, 134, 1, 2, 1, 4, 116, 
	3, 4, 3, 6, 115, 115, 115, 5, 
	6, 5, 8, 97, 7, 8, 7, 10, 
	79, 9, 10, 9, 12, 60, 11, 12, 
	11, 14, 59, 58, 58, 13, 14, 13, 
	16, 40, 15, 16, 15, 18, 21, 21, 
	22, 17, 18, 17, 135, 20, 20, 20, 
	20, 19, 135, 19, 135, 20, 20, 20, 
	20, 19, 18, 22, 17, 18, 23, 17, 
	18, 24, 17, 18, 25, 17, 18, 26, 
	17, 18, 27, 17, 18, 28, 17, 18, 
	29, 17, 18, 30, 17, 18, 31, 17, 
	18, 32, 17, 18, 33, 17, 18, 34, 
	17, 18, 35, 17, 18, 36, 17, 18, 
	37, 17, 18, 38, 17, 18, 39, 17, 
	18, 17, 16, 41, 15, 16, 42, 15, 
	16, 43, 15, 16, 44, 15, 16, 45, 
	15, 16, 46, 15, 16, 47, 15, 16, 
	48, 15, 16, 49, 15, 16, 50, 15, 
	16, 51, 15, 16, 52, 15, 16, 53, 
	15, 16, 54, 15, 16, 55, 15, 16, 
	56, 15, 16, 57, 15, 16, 15, 14, 
	58, 13, 14, 13, 12, 78, 78, 78, 
	78, 78, 61, 78, 78, 11, 12, 78, 
	78, 78, 78, 78, 62, 78, 78, 11, 
	12, 78, 78, 78, 78, 78, 63, 78, 
	78, 11, 12, 78, 78, 78, 78, 78, 
	64, 78, 78, 11, 12, 78, 78, 78, 
	78, 78, 65, 78, 78, 11, 12, 78, 
	78, 78, 78, 78, 66, 78, 78, 11, 
	12, 78, 78, 78, 78, 78, 67, 78, 
	78, 11, 12, 78, 78, 78, 78, 78, 
	68, 78, 78, 11, 12, 78, 78, 78, 
	78, 78, 69, 78, 78, 11, 12, 78, 
	78, 78, 78, 78, 70, 78, 78, 11, 
	12, 78, 78, 78, 78, 78, 71, 78, 
	78, 11, 12, 78, 78, 78, 78, 78, 
	72, 78, 78, 11, 12, 78, 78, 78, 
	78, 78, 73, 78, 78, 11, 12, 78, 
	78, 78, 78, 78, 74, 78, 78, 11, 
	12, 78, 78, 78, 78, 78, 75, 78, 
	78, 11, 12, 78, 78, 78, 78, 78, 
	76, 78, 78, 11, 12, 78, 78, 78, 
	78, 78, 77, 78, 78, 11, 12, 78, 
	78, 78, 78, 78, 78, 78, 11, 12, 
	60, 11, 10, 80, 9, 10, 81, 9, 
	10, 82, 9, 10, 83, 9, 10, 84, 
	9, 10, 85, 9, 10, 86, 9, 10, 
	87, 9, 10, 88, 9, 10, 89, 9, 
	10, 90, 9, 10, 91, 9, 10, 92, 
	9, 10, 93, 9, 10, 94, 9, 10, 
	95, 9, 10, 96, 9, 10, 9, 8, 
	98, 7, 8, 99, 7, 8, 100, 7, 
	8, 101, 7, 8, 102, 7, 8, 103, 
	7, 8, 104, 7, 8, 105, 7, 8, 
	106, 7, 8, 107, 7, 8, 108, 7, 
	8, 109, 7, 8, 110, 7, 8, 111, 
	7, 8, 112, 7, 8, 113, 7, 8, 
	114, 7, 8, 7, 6, 115, 5, 4, 
	117, 3, 4, 118, 3, 4, 119, 3, 
	4, 120, 3, 4, 121, 3, 4, 122, 
	3, 4, 123, 3, 4, 124, 3, 4, 
	125, 3, 4, 126, 3, 4, 127, 3, 
	4, 128, 3, 4, 129, 3, 4, 130, 
	3, 4, 131, 3, 4, 132, 3, 4, 
	133, 3, 4, 3, 2, 134, 134, 1, 
	137, 222, 136, 137, 136, 137, 138, 138, 
	136, 137, 139, 139, 139, 136, 137, 140, 
	136, 137, 141, 144, 181, 184, 187, 201, 
	136, 137, 142, 136, 137, 143, 136, 137, 
	136, 137, 145, 136, 137, 146, 146, 146, 
	146, 167, 146, 146, 136, 137, 147, 136, 
	137, 148, 148, 149, 136, 137, 149, 136, 
	137, 147, 150, 136, 137, 147, 151, 136, 
	137, 147, 152, 136, 137, 147, 153, 136, 
	137, 147, 154, 136, 137, 147, 155, 136, 
	137, 147, 156, 136, 137, 147, 157, 136, 
	137, 147, 158, 136, 137, 147, 159, 136, 
	137, 147, 160, 136, 137, 147, 161, 136, 
	137, 147, 162, 136, 137, 147, 163, 136, 
	137, 147, 164, 136, 137, 147, 165, 136, 
	137, 147, 166, 136, 137, 147, 136, 137, 
	168, 136, 137, 169, 169, 170, 176, 179, 
	175, 136, 137, 170, 176, 175, 136, 137, 
	171, 136, 137, 168, 172, 172, 171, 136, 
	137, 173, 173, 174, 136, 137, 174, 136, 
	137, 168, 174, 136, 137, 168, 170, 172, 
	172, 175, 136, 137, 177, 136, 137, 178, 
	136, 137, 168, 136, 137, 180, 136, 137, 
	178, 136, 137, 182, 136, 137, 183, 183, 
	183, 136, 137, 183, 183, 183, 136, 137, 
	185, 136, 137, 186, 136, 137, 186, 136, 
	137, 188, 136, 137, 189, 189, 190, 196, 
	199, 195, 136, 137, 190, 196, 195, 136, 
	137, 191, 136, 137, 192, 192, 191, 136, 
	137, 193, 193, 194, 136, 137, 194, 136, 
	137, 194, 136, 137, 190, 192, 192, 195, 
	136, 137, 197, 136, 137, 198, 136, 137, 
	136, 137, 200, 136, 137, 198, 136, 137, 
	202, 136, 137, 203, 203, 204, 136, 137, 
	204, 136, 137, 205, 136, 137, 206, 136, 
	137, 207, 136, 137, 208, 136, 137, 209, 
	136, 137, 210, 136, 137, 211, 136, 137, 
	212, 136, 137, 213, 136, 137, 214, 136, 
	137, 215, 136, 137, 216, 136, 137, 217, 
	136, 137, 218, 136, 137, 219, 136, 137, 
	220, 136, 137, 221, 136, 137, 136, 137, 
	222, 136, 0
];

static byte[] _sam_alignment_trans_actions = [
	61, 7, 7, 0, 0, 0, 0, 49, 
	0, 0, 0, 0, 17, 17, 17, 0, 
	0, 0, 0, 49, 0, 0, 0, 0, 
	49, 0, 0, 0, 0, 49, 0, 0, 
	0, 0, 0, 25, 25, 0, 0, 0, 
	0, 49, 0, 0, 0, 0, 1, 1, 
	49, 0, 0, 0, 0, 31, 31, 31, 
	31, 0, 0, 0, 33, 0, 0, 0, 
	0, 0, 0, 49, 0, 52, 3, 0, 
	52, 3, 0, 52, 3, 0, 52, 3, 
	0, 52, 3, 0, 52, 3, 0, 52, 
	3, 0, 52, 3, 0, 52, 3, 0, 
	52, 3, 0, 52, 3, 0, 52, 3, 
	0, 52, 3, 0, 52, 3, 0, 52, 
	3, 0, 52, 3, 0, 52, 3, 0, 
	52, 0, 29, 3, 0, 29, 3, 0, 
	29, 3, 0, 29, 3, 0, 29, 3, 
	0, 29, 3, 0, 29, 3, 0, 29, 
	3, 0, 29, 3, 0, 29, 3, 0, 
	29, 3, 0, 29, 3, 0, 29, 3, 
	0, 29, 3, 0, 29, 3, 0, 29, 
	3, 0, 29, 3, 0, 29, 0, 27, 
	0, 0, 23, 0, 0, 64, 64, 64, 
	64, 64, 3, 64, 64, 0, 0, 64, 
	64, 64, 64, 64, 3, 64, 64, 0, 
	0, 64, 64, 64, 64, 64, 3, 64, 
	64, 0, 0, 64, 64, 64, 64, 64, 
	3, 64, 64, 0, 0, 64, 64, 64, 
	64, 64, 3, 64, 64, 0, 0, 64, 
	64, 64, 64, 64, 3, 64, 64, 0, 
	0, 64, 64, 64, 64, 64, 3, 64, 
	64, 0, 0, 64, 64, 64, 64, 64, 
	3, 64, 64, 0, 0, 64, 64, 64, 
	64, 64, 3, 64, 64, 0, 0, 64, 
	64, 64, 64, 64, 3, 64, 64, 0, 
	0, 64, 64, 64, 64, 64, 3, 64, 
	64, 0, 0, 64, 64, 64, 64, 64, 
	3, 64, 64, 0, 0, 64, 64, 64, 
	64, 64, 3, 64, 64, 0, 0, 64, 
	64, 64, 64, 64, 3, 64, 64, 0, 
	0, 64, 64, 64, 64, 64, 3, 64, 
	64, 0, 0, 64, 64, 64, 64, 64, 
	3, 64, 64, 0, 0, 64, 64, 64, 
	64, 64, 3, 64, 64, 0, 0, 64, 
	64, 64, 64, 64, 64, 64, 0, 21, 
	84, 0, 15, 3, 0, 15, 3, 0, 
	15, 3, 0, 15, 3, 0, 15, 3, 
	0, 15, 3, 0, 15, 3, 0, 15, 
	3, 0, 15, 3, 0, 15, 3, 0, 
	15, 3, 0, 15, 3, 0, 15, 3, 
	0, 15, 3, 0, 15, 3, 0, 15, 
	3, 0, 15, 3, 0, 15, 0, 13, 
	3, 0, 13, 3, 0, 13, 3, 0, 
	13, 3, 0, 13, 3, 0, 13, 3, 
	0, 13, 3, 0, 13, 3, 0, 13, 
	3, 0, 13, 3, 0, 13, 3, 0, 
	13, 3, 0, 13, 3, 0, 13, 3, 
	0, 13, 3, 0, 13, 3, 0, 13, 
	3, 0, 13, 0, 19, 0, 0, 11, 
	3, 0, 11, 3, 0, 11, 3, 0, 
	11, 3, 0, 11, 3, 0, 11, 3, 
	0, 11, 3, 0, 11, 3, 0, 11, 
	3, 0, 11, 3, 0, 11, 3, 0, 
	11, 3, 0, 11, 3, 0, 11, 3, 
	0, 11, 3, 0, 11, 3, 0, 11, 
	3, 0, 11, 0, 9, 0, 0, 0, 
	0, 67, 0, 0, 0, 0, 43, 43, 
	0, 0, 0, 0, 0, 0, 0, 45, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 37, 0, 47, 
	0, 0, 0, 0, 0, 41, 41, 41, 
	41, 41, 41, 41, 0, 0, 0, 0, 
	0, 1, 1, 49, 0, 0, 49, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 3, 0, 
	88, 55, 3, 0, 88, 55, 0, 0, 
	0, 0, 0, 5, 5, 5, 5, 5, 
	5, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 93, 58, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	93, 58, 0, 0, 93, 58, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 93, 58, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 39, 39, 
	39, 0, 73, 0, 0, 0, 0, 0, 
	0, 0, 0, 39, 0, 70, 0, 0, 
	0, 0, 0, 0, 5, 5, 5, 5, 
	5, 5, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 80, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	80, 0, 0, 80, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 80, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 1, 1, 49, 0, 0, 
	49, 0, 76, 3, 0, 76, 3, 0, 
	76, 3, 0, 76, 3, 0, 76, 3, 
	0, 76, 3, 0, 76, 3, 0, 76, 
	3, 0, 76, 3, 0, 76, 3, 0, 
	76, 3, 0, 76, 3, 0, 76, 3, 
	0, 76, 3, 0, 76, 3, 0, 76, 
	3, 0, 76, 3, 0, 76, 0, 0, 
	35, 0, 0
];

static byte[] _sam_alignment_eof_actions = [
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 47, 
	0, 0, 0, 0, 0, 88, 88, 88, 
	88, 88, 88, 88, 88, 88, 88, 88, 
	88, 88, 88, 88, 88, 88, 88, 0, 
	0, 0, 0, 93, 0, 0, 93, 93, 
	0, 0, 93, 0, 0, 0, 0, 73, 
	0, 0, 70, 0, 0, 0, 0, 80, 
	0, 0, 80, 80, 0, 0, 80, 0, 
	0, 0, 0, 0, 76, 76, 76, 76, 
	76, 76, 76, 76, 76, 76, 76, 76, 
	76, 76, 76, 76, 76, 76, 0
];

static int sam_alignment_start = 0;
static int sam_alignment_first_final = 135;
static int sam_alignment_error = -1;

static int sam_alignment_en_alignment = 0;


#line 234 "sam_alignment.rl"


import bio.sam.header;
import bio.bam.read;
import bio.bam.tagvalue;
import bio.bam.utils.tagstoragebuilder;

import std.array;
import std.conv;
import std.typecons;
import std.outbuffer;
import std.c.stdlib;

class AlignmentBuildStorage {
    Appender!(CigarOperation[]) cigar_appender;
    OutBuffer outbuffer;
    TagStorageBuilder tag_storage_builder;

    this() {
        cigar_appender = appender!(CigarOperation[])();
        outbuffer = new OutBuffer();
        tag_storage_builder = TagStorageBuilder.create();
    }

    void clear() {
        cigar_appender.clear();
        tag_storage_builder.clear();
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
    }
}

BamRead parseAlignmentLine(string line, SamHeader header, 
                             AlignmentBuildStorage b=null) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    if (b is null) {
        b = new AlignmentBuildStorage();
    } else {
        b.clear();
    }

    byte current_sign = 1;

    size_t read_name_beg; // position of beginning of QNAME
    size_t read_name_end; // position past the end of QNAME

    size_t sequence_beg; // position of SEQ start
    string sequence;     // SEQ

    uint cigar_op_len;   // length of CIGAR operation
    char cigar_op_chr;   // CIGAR operation

    size_t cigar_op_len_start; // position of start of CIGAR operation
    
    auto cigar = b.cigar_appender;

    long int_value;                      // for storing temporary integers
    float float_value;                   // for storing temporary floats
    size_t float_beg;                    // position of start of current float
    auto outbuffer = b.outbuffer;        // used to build tag values which hold arrays
    char arraytype;                      // type of last array tag value

    ushort flag;
    uint pos;
    uint mate_pos;
    ubyte mapping_quality; 
    int template_length;
    ubyte* qual_ptr = null;
    size_t qual_index; 

    string current_tag;
    Value current_tagvalue;

    size_t tag_key_beg, tagvalue_beg;
    size_t rname_beg, rnext_beg;

    int ref_id = -1;
    int mate_ref_id = -1;
    
    auto builder = b.tag_storage_builder;

    
#line 624 "sam_alignment.d"
	{
	cs = sam_alignment_start;
	}

#line 320 "sam_alignment.rl"
    
#line 631 "sam_alignment.d"
	{
	int _klen;
	uint _trans;
	byte* _acts;
	uint _nacts;
	char* _keys;

	if ( p == pe )
		goto _test_eof;
_resume:
	_keys = &_sam_alignment_trans_keys[_sam_alignment_key_offsets[cs]];
	_trans = _sam_alignment_index_offsets[cs];

	_klen = _sam_alignment_single_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += cast(uint)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _sam_alignment_range_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += cast(uint)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	cs = _sam_alignment_trans_targs[_trans];

	if ( _sam_alignment_trans_actions[_trans] == 0 )
		goto _again;

	_acts = &_sam_alignment_actions[_sam_alignment_trans_actions[_trans]];
	_nacts = cast(uint) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 23 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	break;
	case 1:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
	break;
	case 2:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	break;
	case 3:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
	break;
	case 4:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	break;
	case 5:
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
	break;
	case 6:
#line 43 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
	break;
	case 7:
#line 44 "sam_alignment.rl"
	{ read_name_end = p - line.ptr; }
	break;
	case 8:
#line 48 "sam_alignment.rl"
	{ flag = to!ushort(int_value); }
	break;
	case 9:
#line 49 "sam_alignment.rl"
	{ pos = to!uint(int_value); }
	break;
	case 10:
#line 50 "sam_alignment.rl"
	{ mapping_quality = to!ubyte(int_value); }
	break;
	case 11:
#line 54 "sam_alignment.rl"
	{ rname_beg = p - line.ptr; }
	break;
	case 12:
#line 55 "sam_alignment.rl"
	{
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]); 
    }
	break;
	case 13:
#line 64 "sam_alignment.rl"
	{ cigar_op_len = to!uint(int_value); }
	break;
	case 14:
#line 65 "sam_alignment.rl"
	{ cigar_op_chr = (*p); }
	break;
	case 15:
#line 66 "sam_alignment.rl"
	{ 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }
	break;
	case 16:
#line 73 "sam_alignment.rl"
	{
        mate_ref_id = ref_id;
    }
	break;
	case 17:
#line 77 "sam_alignment.rl"
	{ rnext_beg = p - line.ptr; }
	break;
	case 18:
#line 78 "sam_alignment.rl"
	{
        mate_ref_id = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
    }
	break;
	case 19:
#line 85 "sam_alignment.rl"
	{ mate_pos = to!uint(int_value); }
	break;
	case 20:
#line 86 "sam_alignment.rl"
	{ template_length = to!int(int_value); }
	break;
	case 21:
#line 91 "sam_alignment.rl"
	{ sequence_beg = p - line.ptr; }
	break;
	case 22:
#line 92 "sam_alignment.rl"
	{ sequence = line[sequence_beg .. p - line.ptr]; }
	break;
	case 23:
#line 97 "sam_alignment.rl"
	{
        if (sequence.length > 1024) {
            qual_ptr = (new ubyte[sequence.length]).ptr;
        } else {
            qual_ptr = cast(ubyte*)alloca(sequence.length);
            if (!qual_ptr) {
                qual_ptr = (new ubyte[sequence.length]).ptr;
            }
        }
        qual_index = 0;
    }
	break;
	case 24:
#line 109 "sam_alignment.rl"
	{
        qual_ptr[qual_index++] = cast(ubyte)((*p) - 33);
    }
	break;
	case 25:
#line 126 "sam_alignment.rl"
	{ current_tagvalue = Value((*p)); }
	break;
	case 26:
#line 127 "sam_alignment.rl"
	{ 
        if (int_value < 0) {
            if (int_value >= byte.min) {
                current_tagvalue = Value(to!byte(int_value));
            } else if (int_value >= short.min) {
                current_tagvalue = Value(to!short(int_value));
            } else if (int_value >= int.min) {
                current_tagvalue = Value(to!int(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                current_tagvalue = Value(to!ubyte(int_value));
            } else if (int_value <= ushort.max) {
                current_tagvalue = Value(to!ushort(int_value));
            } else if (int_value <= uint.max) {
                current_tagvalue = Value(to!uint(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
	break;
	case 27:
#line 151 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	break;
	case 28:
#line 153 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(float_value);
    }
	break;
	case 29:
#line 157 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]); 
    }
	break;
	case 30:
#line 161 "sam_alignment.rl"
	{
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]);
        current_tagvalue.setHexadecimalFlag();
    }
	break;
	case 31:
#line 170 "sam_alignment.rl"
	{
        // it might be not the best idea to use outbuffer;
        // the better idea might be two-pass approach
        // when first pass is for counting commas, and
        // the second is for filling allocated array
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
        arraytype = (*p);
    }
	break;
	case 32:
#line 180 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }
	break;
	case 33:
#line 193 "sam_alignment.rl"
	{ 
        outbuffer.write(float_value); 
    }
	break;
	case 34:
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
	break;
	case 35:
#line 223 "sam_alignment.rl"
	{ tag_key_beg = p - line.ptr; }
	break;
	case 36:
#line 224 "sam_alignment.rl"
	{ current_tag = line[tag_key_beg .. p - line.ptr]; }
	break;
	case 37:
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
#line 937 "sam_alignment.d"
		default: break;
		}
	}

_again:
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	byte* __acts = &_sam_alignment_actions[_sam_alignment_eof_actions[cs]];
	uint __nacts = cast(uint) *__acts++;
	while ( __nacts-- > 0 ) {
		switch ( *__acts++ ) {
	case 3:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
	break;
	case 5:
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
	break;
	case 26:
#line 127 "sam_alignment.rl"
	{ 
        if (int_value < 0) {
            if (int_value >= byte.min) {
                current_tagvalue = Value(to!byte(int_value));
            } else if (int_value >= short.min) {
                current_tagvalue = Value(to!short(int_value));
            } else if (int_value >= int.min) {
                current_tagvalue = Value(to!int(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                current_tagvalue = Value(to!ubyte(int_value));
            } else if (int_value <= ushort.max) {
                current_tagvalue = Value(to!ushort(int_value));
            } else if (int_value <= uint.max) {
                current_tagvalue = Value(to!uint(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
	break;
	case 28:
#line 153 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(float_value);
    }
	break;
	case 29:
#line 157 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]); 
    }
	break;
	case 30:
#line 161 "sam_alignment.rl"
	{
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]);
        current_tagvalue.setHexadecimalFlag();
    }
	break;
	case 32:
#line 180 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }
	break;
	case 33:
#line 193 "sam_alignment.rl"
	{ 
        outbuffer.write(float_value); 
    }
	break;
	case 34:
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
	break;
	case 37:
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
#line 1047 "sam_alignment.d"
		default: break;
		}
	}
	}

	}

#line 321 "sam_alignment.rl"

    auto read = BamRead(line[read_name_beg .. read_name_end], 
                        sequence,
                        cigar.data,
                        builder.data);

    if (qual_ptr !is null && qual_index == sequence.length) {
        read.base_qualities = qual_ptr[0 .. sequence.length];
    }

    read.flag = flag;
    read.mapping_quality = mapping_quality;
    read.position = pos - 1; // we use 0-based coordinates, not 1-based
    read.template_length = template_length;
    read.mate_position = mate_pos - 1; // also 0-based
    read.ref_id = ref_id;
    read.mate_ref_id = mate_ref_id;

    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto header = new SamHeader("@SQ\tSN:20\tLN:1234567");
    auto alignment = parseAlignmentLine(line, header);
    assert(alignment.name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60032);
    assert(alignment.mapping_quality == 25);
    assert(alignment.mate_position == 60032);
    assert(alignment.ref_id == 0);
    assert(alignment.mate_ref_id == 0);
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
    assert(to!char(alignment["XT"]) == 'U');

    import std.stdio;
    import bio.bam.serialization.sam;
    import bio.bam.reference;

    auto info = ReferenceSequenceInfo("20", 1234567);

    auto invalid_cigar_string = "1\t100\t20\t50000\t30\tMZABC\t=\t50000\t0\tACGT\t####";
    alignment = parseAlignmentLine(invalid_cigar_string, header);
    assert(equal(alignment.sequence(), "ACGT"));

    auto invalid_tag_and_qual = "2\t100\t20\t5\t40\t27M30X5D\t=\t3\t10\tACT\t !\n\tX1:i:7\tX3:i:zzz\tX4:i:5";
    alignment = parseAlignmentLine(invalid_tag_and_qual, header);
    assert(alignment.base_qualities == [255, 255, 255]); // i.e. invalid
    assert(to!ubyte(alignment["X1"]) == 7);
    assert(alignment["X3"].is_nothing);
    assert(to!ubyte(alignment["X4"]) == 5);

}
