# pyaces DONE

완료한 작업 기록. 새 기록은 맨 위에 추가한다(최신순).
[TODO.md](TODO.md)에서 온 항목은 같은 번호(`T-<번호>`)를 쓰고, TODO에 없던 작업은 번호 없이 `D-<날짜>`로 적는다.

## 작성 규칙

```markdown
### T-00. 제목 (YYYY-MM-DD)
- **커밋:** `해시` (여러 개면 모두)
- **한 일:** 무엇을 어떻게 바꿨나
- **검증:** 로컬 / ccsd 결과 (aces2, pyaces, 링크 감사), 추가 확인
- **남은 것 / 메모:** 후속 작업이나 주의점 (없으면 생략)
```

---

### T-16. 테스트 요약 누락과 분류기 문제 (2026-10-23)
- **커밋:** `95cc01a2`
- **한 일:**
  - `run_tests.sh`의 pyaces 요약이 러너 Summary 블록의 모든 `STATUS: n` 줄을 읽는다. PASS → PASS, MISMATCH·FAIL_CLEAN → FAIL, 나머지(CRASH, NO_RESULT, TIMEOUT …) → OTHER.
  - TSV에 손으로 고쳤던 규칙을 `classify()`로 옮겼다: `SCF_TYPE=KS` → KS/Runintgrt, `FOCK=AO` MBPT(2) → Rundirmp2, CCSD(TQ)·`BRUECKNER`·`HFSTABILITY=FOLLOW` 제외. `classify_suite.py --check`(TSV가 규칙과 다르면 exit 1)를 추가했다.
  - TSV를 재생성했다. 실행 대상 51행은 그대로이고, 제외 행의 설명 문구 3개와 020의 iuhf(REF=UHF라 0 → 1)만 바뀌었다.
  - 러너의 기본 실행 디렉터리를 저장소의 `build_logs/pyaces_suite_run`으로 바꾸고, HiPerGator에서는 `run_tests.sh`가 `/blue` 경로를 지정하게 했다.
- **검증:**
  - 가짜 로그로 비교: MISMATCH 2·TIMEOUT 1·CRASH 1 → 새 요약 `TOTAL=4 FAIL=2 OTHER=2`(이전 코드는 모두 `?`). PASS 없이 CRASH 3 → `TOTAL=3 OTHER=3`(이전 `?`).
  - `classify_suite.py --check`: 로컬·ccsd 모두 일치(125행).
  - 로컬·ccsd 스위트 결과는 T-03 항목 참고.
- **메모:** 시험 중 잘못 짠 명령이 로컬 `build_logs/unified_test_pyaces.log`를 덮어써서 pyaces 테스트를 다시 돌려 복구했다(51/51).

### T-07. Python 버전 고정 해제 (2026-10-23)
- **커밋:** `95cc01a2`
- **한 일:** `build.sh`가 `sysconfig.get_config_var('EXT_SUFFIX')`로 확장 모듈 이름을 정한다(하드코딩된 `cpython-312` 제거). 감사 도구의 `.so` 탐색은 원래 `aces2py*.so`라 수정이 필요 없었다.
- **검증:** 3.12에서 이름(`aces2py.cpython-312-x86_64-linux-gnu.so`)과 결과 동일(로컬·ccsd).
- **남은 것:** 다른 Python 버전에서는 시험하지 못했다. 로컬에는 3.12만 있고, ccsd의 python/3.10·3.11·3.6에는 numpy가 없다.

### T-04. 이름 바꾸기 실패 시 빌드 중단 (2026-10-23)
- **커밋:** `95cc01a2`
- **한 일:** 패치 전에 `EXTRA_LOCALIZE_SYMBOLS`의 모든 쌍(자동 추가분 포함)을 검사한다. 모듈이 `MAIN_COLLISION_MODULES`/`EXTRA_LOCALIZE_MODULES`에 없거나 아카이브에 `T <symbol>`이 없으면 목록을 출력하고 빌드를 멈춘다. `objcopy --redefine-sym`의 `2>/dev/null || true`도 없앴다.
- **검증:**
  - 현재 60쌍은 모두 유효(사전 점검).
  - 일부러 `vcc:nosuch_`, `notamod:finish_`를 넣은 빌드: 두 쌍을 모두 보고하고 rc=1로 멈춤.

### T-03. COMMON 블록에 연결된 호출 검사 (2026-10-23)
- **커밋:** `95cc01a2`
- **한 일:** `collision_audit.py --calls`를 추가했다. pyaces Fortran 오브젝트(`Main/build/*.o`, C 오브젝트 2개 제외)의 `U` 심볼이 `.so`에서 데이터(B/b/D/d)이면 보고하고 exit 1. 빌드 5단계(import 시험 직후)에서 실행하며, 걸리면 빌드가 실패한다.
- **검증:**
  - f97a9759 이전의 `runhbar.F`/`buildhbar.F`로 오브젝트를 만들어 검사: `hbar_` 2건 보고, rc=1. 현재 오브젝트: 0건, rc=0.
  - 로컬: pyaces 빌드 OK, aces2 125/125, pyaces 51/51, 감사 296(기준선 일치), `--calls` 0건.
  - ccsd (PBS 7713/7714): 같은 결과.

### D-2026-10-23b. pyaces 문서 작성 (2026-10-23)
- **커밋:** `86a51cda`, `3a123417`, `daecf899`, `259a4fda`
- **한 일:**
  - 개발자 매뉴얼을 작성했다(`docs/developer_manual.html`, `docs/pyaces_developer_manual.pdf`, 25쪽, 영어).
  - 개선 과제 문서를 한국어·영어로 작성했다(`docs/pyaces_improvements_ko.pdf`, `docs/pyaces_improvements_en.pdf`).
  - 개선 과제를 [TODO.md](TODO.md)로 정리하고 이 완료 기록을 시작했다.
- **검증:** PDF 생성·렌더링 확인, 영어판에 한글이 남지 않았음을 확인.

### D-2026-10-23a. 매뉴얼 작성 중 찾은 결함 3개 수정 (2026-10-23)
- **커밋:** `f97a9759`
- **한 일:**
  - `Runhbar`/`Buildhbar`의 `Call Hbar` → `Call Hbar_`. `hbar_`는 서브루틴이 아니라 lambda의 `COMMON /HBAR/`(.so에서 `B hbar_`)이어서, SCF·CC가 끝난 상태의 분기가 BSS로 점프할 수 있었다.
  - `run*.F`의 `Integer Iuhf` → `Integer*8 Iuhf`. f2py가 4바이트로 넘기고 `-i8` Fortran이 8바이트를 읽던 불일치를 없앴다.
  - `Getreclenpy`, `Cgetrecpy`에 `Cf2py intent(out)`을 붙여 길이와 문자열을 Python으로 반환하게 했다. `Cgetrecpy`의 디버그 출력은 제거했다.
- **검증:**
  - 로컬: pyaces 51/51, 감사 296(기준선 일치). 케이스 026에서 `Init → Runscf → Runhbar` 실행, CCSD 에너지가 기준값과 7×10⁻¹¹ 이내.
  - ccsd (PBS 7700/7701): pyaces 51/51, 감사 296.
- **API 변경:** `Length = Getreclenpy(Name)`, `Dest = Cgetrecpy(Name, Length)`.
- **메모:** `Runhbar`를 `Runcc` 다음에 부르면 SCF 분기로 들어가 Vtran을 두 번 부른다. 매뉴얼에 명시했다.

### D-2026-10-14b. 링크 결과가 빌드마다 바뀌던 문제 수정 (2026-10-14)
- **커밋:** `eab3c4e6`
- **한 일:**
  - aces2 makefile의 소스 목록을 `$(sort …)`로 정렬했다(GNU Make 4.2는 wildcard를 정렬하지 않고, 4.3은 로케일 순서로 정렬).
  - `pyaces/build.sh`에 `export LC_ALL=C`를 넣었다.
  - 링크 순서에 따라 승자가 바뀌던 `getdt1`(vcc), `fndmaxs`(vee), `t2l2rng`·`r1l2y2c`(lcct)를 이름 바꾸기 목록에 추가했다. lcct가 vee의 `r1l2y2c`를 쓰던 기준선 항목이 해소되어 기준선이 297 → 296이 됐다.
- **검증:** 로컬·ccsd 모두 aces2 125/125, pyaces 51/51, 감사 296으로 일치.

### D-2026-10-14a. ccsd 서버 설치와 oneAPI 2023.2 대응 (2026-10-14)
- **커밋:** `73944a6d`, `dbed9377`
- **한 일:**
  - ccsd(`/home/ycpark/program/ACES`)에 설치. 갱신은 `git push ccsd master`, 빌드·테스트는 PBS로만.
  - `pyaces/Main/init.F`의 `#include`를 1열로 옮겼다(ifort 2023.2가 들여쓴 include를 거부).
  - aces2 Intel 빌드에 `-assume dummy_aliases`를 추가했다(ifx 2023.2가 `-O2`에서 vmol의 `VGAMMA2`를 잘못 컴파일해 SCF 에너지가 0.2 Eh 틀렸음).
- **검증:** ccsd에서 aces2 125/125, aces3 37/37, pyaces 51/51.
