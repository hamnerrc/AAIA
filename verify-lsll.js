#!/usr/bin/env node
/**
 * Verify that reco_lsll.txt segments are correct:
 * - Scramble + inspection + path gets to a 3-pair state
 * - LSLL algorithm completes the 4th pair
 */
'use strict';
const fs = require('fs');
const path = require('path');
const { parseProReferences } = require('./pro-references.js');
const { applyAlgorithm, SOLVED_FACELETS, canonicalizeForEngine } = require('../js/facelet-cube.js');
const { solvedFlags } = require('../js/facelet-flags.js');

function getPairCount(facelets) {
  const f = solvedFlags(facelets);
  const slots = ['BL', 'BR', 'FL', 'FR'];
  return slots.filter(s => f[s]).length;
}

function verifyLSLL(file) {
  const solves = parseProReferences(fs.readFileSync(file, 'utf8'));
  const results = { valid: 0, invalid: 0, errors: [] };

  for (let i = 0; i < solves.length; i++) {
    const solve = solves[i];
    try {
      // Apply scramble and inspection
      let facelets = applyAlgorithm(SOLVED_FACELETS, [solve.scramble, solve.inspection].filter(Boolean).join(' '));

      // In reco_lsll format: step 0 is path to 3-pair, step 1 is LSLL
      if (solve.steps.length < 2) {
        results.errors.push(`[${i}] ${solve.source}: expected 2 steps, got ${solve.steps.length}`);
        results.invalid++;
        continue;
      }

      const pathStep = solve.steps[0];
      const lsllStep = solve.steps[1];

      // Apply path to 3-pair state
      facelets = applyAlgorithm(facelets, pathStep.alg);
      const pairsAfterPath = getPairCount(facelets);

      if (pairsAfterPath !== 3) {
        results.errors.push(`[${i}] ${solve.source}: path ends with ${pairsAfterPath} pairs, expected 3`);
        results.invalid++;
        continue;
      }

      // Apply LSLL algorithm
      facelets = applyAlgorithm(facelets, lsllStep.alg);
      const pairsAfterLSLL = getPairCount(facelets);

      if (pairsAfterLSLL !== 4) {
        results.errors.push(`[${i}] ${solve.source}: LSLL ends with ${pairsAfterLSLL} pairs, expected 4`);
        results.invalid++;
        continue;
      }

      results.valid++;
    } catch (e) {
      results.errors.push(`[${i}] ${solve.source}: ${e.message}`);
      results.invalid++;
    }
  }

  console.log(`Verified ${results.valid} valid LSLL segments, ${results.invalid} invalid`);
  if (results.errors.length > 0 && results.errors.length <= 20) {
    console.log('Errors:');
    results.errors.forEach(e => console.log('  ' + e));
  } else if (results.errors.length > 20) {
    console.log(`First 20 errors:`);
    results.errors.slice(0, 20).forEach(e => console.log('  ' + e));
  }

  return results.invalid === 0 ? 0 : 1;
}

if (require.main === module) {
  const file = process.argv[2] || path.join(__dirname, '..', 'data', 'reco_lsll.txt');
  process.exit(verifyLSLL(file));
}

module.exports = { verifyLSLL };
